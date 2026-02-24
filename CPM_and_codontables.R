# ============================
# CODON COUNTS and CPM tables
# ============================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(tibble)
})

# ============================
# 1) Compute CPM for one file
# ============================

counts_path <- "/Users/nross1994/Documents/Doore Lab/RNAseq Analysis/examples for nykki/LB_37_A1-A11/A9.txt"
count_col   <- "A9.bam"

counts <- readr::read_tsv(counts_path, show_col_types = FALSE)

stopifnot("Geneid" %in% names(counts))
stopifnot(count_col %in% names(counts))

total_counts <- sum(counts[[count_col]], na.rm = TRUE)

# Build CPM table keyed by a JOIN ID
# - host: S_RS....
# - phage in counts: FDI99_gp### (or similar) often prefixed "FDI99_"
expr <- counts %>%
  transmute(
    geneid_raw     = .data$Geneid,
    locus_tag_raw  = sub("^gene-", "", .data$Geneid),
    CPM            = .data[[count_col]] / total_counts * 1e6
  ) %>%
  mutate(
    locus_tag_join = sub("^FDI99_", "", locus_tag_raw)  # normalize for joining
  ) %>%
  group_by(locus_tag_join) %>%
  summarise(CPM = sum(CPM, na.rm = TRUE), .groups = "drop")

cat("Counts rows:", nrow(counts), "\n")
cat("Expr distinct locus_tag_join:", n_distinct(expr$locus_tag_join), "\n")


# =======================================
# 2) Parse GenBank and extract CDS sequences
# =======================================

host_path <- "/Users/nross1994/Documents/Doore Lab/RNAseq Analysis/examples for nykki/2457T_genbank.gb"
phage_path <- "/Users/nross1994/Documents/Doore Lab/RNAseq Analysis/examples for nykki/Sf14_genbank.gb"

host_txt <- readLines(host_path)
phage_txt <- readLines(phage_path)

extract_genome_from_gb <- function(txt) {
  origin_i <- grep("^ORIGIN", txt)[1]
  end_i <- grep("^//", txt)
  end_i <- end_i[end_i > origin_i][1]
  
  seq_lines <- txt[(origin_i + 1):(end_i - 1)]
  genome <- paste(seq_lines, collapse = "")
  genome <- gsub("[^acgtACGT]", "", genome)
  tolower(genome)
}

host_genome  <- extract_genome_from_gb(host_txt)
phage_genome <- extract_genome_from_gb(phage_txt)

cat("Host genome length:", nchar(host_genome), "\n")
cat("Phage genome length:", nchar(phage_genome), "\n")

get_feat_block <- function(txt) {
  feat_i   <- which(txt == "FEATURES             Location/Qualifiers")[1]
  origin_i <- grep("^ORIGIN", txt)[1]
  txt[(feat_i + 1):(origin_i - 1)]
}

is_feature_start <- function(line) {
  grepl("^\\s{5}\\S", line) && !grepl("^\\s{5}/", line)
}

split_records <- function(feature_lines) {
  starts <- which(vapply(feature_lines, is_feature_start, logical(1)))
  ends   <- c(starts[-1] - 1, length(feature_lines))
  Map(function(s, e) feature_lines[s:e], starts, ends)
}

parse_feature_header <- function(line) {
  line <- sub("^\\s+", "", line)
  key <- sub("\\s+.*$", "", line)
  loc <- sub("^\\S+\\s+", "", line)
  list(type = key, location = loc)
}

get_qual <- function(rec_text, key) {
  m1 <- stringr::str_match(rec_text, paste0("/", key, '="([^"]*)"'))
  if (!is.na(m1[,2])) return(m1[,2])
  m2 <- stringr::str_match(rec_text, paste0("/", key, "=([^\\s]+)"))
  if (!is.na(m2[,2])) return(m2[,2])
  NA_character_
}

get_qual_num <- function(rec_text, key) {
  x <- get_qual(rec_text, key)
  suppressWarnings(as.integer(x))
}

# build feature table
build_feature_df <- function(records) {
  purrr::map_df(records, function(rec) {
    hdr <- parse_feature_header(rec[1])
    rec_text <- paste(rec, collapse = "\n")
    tibble::tibble(
      type         = hdr$type,
      location     = hdr$location,
      gene         = get_qual(rec_text, "gene"),
      locus_tag    = get_qual(rec_text, "locus_tag"),
      product      = get_qual(rec_text, "product"),
      protein_id   = get_qual(rec_text, "protein_id"),
      codon_start  = get_qual_num(rec_text, "codon_start"),
      transl_table = get_qual_num(rec_text, "transl_table"),
      note         = get_qual(rec_text, "note")
    )
  })
}

# Build features
host_records  <- split_records(get_feat_block(host_txt))
phage_records <- split_records(get_feat_block(phage_txt))

host_features  <- build_feature_df(host_records)
phage_features <- build_feature_df(phage_records)

cat("Host feature types:\n")
print(dplyr::count(host_features, type) %>% dplyr::arrange(dplyr::desc(n)), n = 50)

cat("Phage feature types:\n")
print(dplyr::count(phage_features, type) %>% dplyr::arrange(dplyr::desc(n)), n = 50)

# ----------------------------
# Sequence extraction helpers
# ----------------------------
rev_comp_str <- function(s) {
  s <- tolower(s)
  s <- chartr("acgtn", "tgcan", s)
  paste(rev(strsplit(s, "")[[1]]), collapse = "")
}

extract_by_ranges <- function(genome, ranges) {
  parts <- vapply(ranges, function(r) substr(genome, r[1], r[2]), character(1))
  paste(parts, collapse = "")
}

parse_location_to_ranges <- function(loc) {
  loc <- gsub("\\s", "", loc)
  m <- stringr::str_match_all(loc, "([<>]?\\d+)\\.\\.([<>]?\\d+)")[[1]]
  if (nrow(m) == 0) stop("No ranges found in location: ", loc)
  to_int <- function(x) as.integer(gsub("[<>]", "", x))
  ranges <- lapply(seq_len(nrow(m)), function(k) c(to_int(m[k,2]), to_int(m[k,3])))
  list(
    ranges = ranges,
    is_complement = grepl("^complement\\(", loc)
  )
}

extract_feature_seq <- function(genome, loc) {
  parsed <- parse_location_to_ranges(loc)
  s <- extract_by_ranges(genome, parsed$ranges)
  if (parsed$is_complement) s <- rev_comp_str(s)
  s
}

# Add sequences
host_features <- host_features %>%
  mutate(seq = purrr::map_chr(location, ~extract_feature_seq(host_genome, .x)))

phage_features <- phage_features %>%
  mutate(seq = purrr::map_chr(location, ~extract_feature_seq(phage_genome, .x)))


# ======================================
# 3) Build CDS table (KEY STEP FOR 4495)
# ======================================

host_cds <- host_features %>%
  dplyr::filter(type == "CDS") %>%
  dplyr::select(type, location, gene, locus_tag, protein_id, product, codon_start, transl_table, note, seq) %>%
  dplyr::mutate(genome = "2457T")

phage_cds <- phage_features %>%
  dplyr::filter(type == "CDS") %>%
  dplyr::filter(!is.na(locus_tag) & locus_tag != "") %>%
  dplyr::mutate(genome = "Sf14")

cat("Host CDS:", nrow(host_cds), "\n")
cat("Phage CDS:", nrow(phage_cds), "\n")

cds_tbl <- dplyr::bind_rows(host_cds, phage_cds) %>%
  dplyr::group_by(genome) %>%
  dplyr::mutate(
    cds_uid     = paste0(genome, "_CDS_", dplyr::row_number()),
    codon_start = dplyr::if_else(is.na(codon_start), 1L, codon_start),
    seq_cds     = substr(seq, codon_start, nchar(seq)),
    locus_tag_join = dplyr::if_else(genome == "Sf14", sub("^FDI99_", "", locus_tag), locus_tag),
    
    # NEW: keep protein_id (host from GenBank; phage will be NA)
    protein_id = dplyr::if_else(genome == "Sf14", NA_character_, protein_id)
  ) %>%
  dplyr::ungroup()

print(dplyr::count(cds_tbl, genome), n = 10)


# ======================================
# 4) Codon counts per CDS
# ======================================

split_codons <- function(seq) {
  if (is.na(seq)) return(character(0))
  seq <- gsub("\\s+", "", seq)
  seq <- tolower(seq)
  L <- nchar(seq)
  if (L < 3) return(character(0))
  seq <- substr(seq, 1, L - (L %% 3))
  substring(seq, seq(1, nchar(seq), by = 3),
            seq(3, nchar(seq), by = 3))
}

count_codons <- function(seq) {
  codons <- split_codons(seq)
  if (length(codons) == 0) {
    return(tibble(codon = character(0), count = integer(0)))
  }
  tab <- table(codons)
  tibble(codon = names(tab), count = as.integer(tab))
}

codon_counts_long <- cds_tbl %>%
  mutate(codon_tbl = purrr::map(seq_cds, count_codons)) %>%
  tidyr::unnest(codon_tbl) %>%
  transmute(
    genome, cds_uid, locus_tag, protein_id, locus_tag_join, gene, product,
    codon = tolower(gsub("\\s+", "", codon)),
    count
  ) %>%
  filter(nchar(codon) == 3) %>%
  filter(str_detect(codon, "^[acgt]{3}$")) %>%
  filter(!codon %in% c("taa","tag","tga"))

cat("Distinct CDS in codon_counts_long:", n_distinct(codon_counts_long$cds_uid), "\n")


# ======================================
# 5) Join CPM weights (missing genes => CPM=0)
# ======================================

codon_counts_weighted <- codon_counts_long %>%
  left_join(expr, by = "locus_tag_join") %>%
  mutate(CPM = coalesce(CPM, 0))

cat("Distinct CDS in codon_counts_weighted:", n_distinct(codon_counts_weighted$cds_uid), "\n")
cat("Genes with CPM>0:", codon_counts_weighted %>% distinct(cds_uid, CPM) %>% summarise(sum(CPM>0)) %>% pull(), "\n")

write.csv(
  codon_counts_weighted,
  "/Users/nross1994/Documents/Doore Lab/RNAseq Analysis/examples for nykki/all_codon_table.csv",
  row.names = FALSE
)
