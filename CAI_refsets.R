# ============================================================
# Build CAI reference sets from merged_counts CSV (host + phage)
#   - Host: 50 genes across ALL samples
#   - Phage: 20 genes across INFECTED samples only
# ============================================================
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.22")

BiocManager::install("edgeR")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(edgeR)
})

# ----------------------------
# LOAD COUNTS AND DESIGN
# ----------------------------
merged_counts_csv <- "/Users/nross1994/Documents/Doore Lab/RNAseq Analysis/merged_counts_ALL.csv"
design_tsv <- "/Users/nross1994/Documents/Doore Lab/RNAseq Analysis/design.tsv"

# if CSV is "wide format" - one column for geneid and one column for each sample - convert to long
wide <- read_csv(merged_counts_csv, show_col_types = FALSE)

# Pivot: all columns except GeneID become sample/count pairs
merged_long <- wide %>%
  pivot_longer(
    cols      = -GeneID,
    names_to  = "sample",
    values_to = "count"
  ) %>%
  mutate(
    GeneID = as.character(GeneID),
    sample = as.character(sample),
    count  = as.numeric(count),
    genome = case_when(
      str_detect(GeneID, "FDI99") ~ "Sf14",     # phage
      TRUE                       ~ "2457T"     # host (fallback)
    )
  )

gene_col   <- "Geneid"     # e.g., "locus_tag", "gene_id", "feature_id"
sample_col <- "sample"        # e.g., "sample", "well", "sample_id"
count_col  <- "count"         # e.g., "count", "counts", "raw_count"
genome_col <- "genome"        # e.g., "genome", "organism", "species"

host_genome_value  <- "2457T" 
phage_genome_value <- "Sf14"  

# read design and isolate infected samples for phage analysis
design <- read_csv(design_tsv, show_col_types = FALSE)
head(design)

design <- mutate(design, sample = as.character(sample), condition = as.character(condition))

infected_samples <- design %>%
  filter(condition == "I") %>%
  pull(sample)

if (length(infected_samples) == 0) stop("No infected samples found (condition == 'I').")

# ---------------
# SET PARAMETERS
# ---------------

# Reference set sizes
n_host  <- 50
n_phage <- 20

# Expression / prevalence thresholds
host_present_cpm <- 1     # gene considered "present" in a sample if CPM >= this
phage_present_cpm <- 1
host_prevalence_min   <- 0.85  # host gene must be present in >=80% samples
phage_prevalence_min  <- 0.70  # phage gene must be present in >=70% INFECTED samples

# Variability filtering (keep "stable" genes):
# we keep genes with sd_logCPM <= sd_quantile_cut among candidates (lower = stricter)
sd_quantile_cut <- 0.50        # 0.50 = keep lower 50% SD (stable half)

# Scoring: higher median is good; higher SD is bad
sd_penalty <- 0.5              # increase to penalize variability more strongly

# Output folder
out_dir <- "./CAI_reference_sets"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Helpers
# ----------------------------
make_logcpm_matrix <- function(df_long, gene_col, sample_col, count_col) {
  mat_df <- df_long %>%
    select(all_of(c(gene_col, sample_col, count_col))) %>%
    mutate(!!count_col := as.numeric(.data[[count_col]])) %>%
    group_by(.data[[gene_col]], .data[[sample_col]]) %>%
    summarise(count = sum(.data[[count_col]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = all_of(sample_col), values_from = count, values_fill = 0)
  
  gene_ids <- mat_df[[gene_col]]
  mat <- as.matrix(mat_df[, setdiff(names(mat_df), gene_col)])
  rownames(mat) <- gene_ids
  
  logcpm <- edgeR::cpm(mat, log = TRUE, prior.count = 1)
  cpm    <- edgeR::cpm(mat, log = FALSE)
  
  list(raw = mat, cpm = cpm, logcpm = logcpm)
}

summarise_gene_stability <- function(cpm_mat, logcpm_mat, present_cpm = 1, sd_penalty = 0.5) {
  present <- cpm_mat >= present_cpm
  
  tibble(
    gene = rownames(logcpm_mat),
    median_logCPM = apply(logcpm_mat, 1, median, na.rm = TRUE),
    sd_logCPM     = apply(logcpm_mat, 1, sd, na.rm = TRUE),
    iqr_logCPM    = apply(logcpm_mat, 1, IQR, na.rm = TRUE),
    prevalence    = rowMeans(present, na.rm = TRUE),
    n_samples     = ncol(logcpm_mat)
  ) %>%
    mutate(score = median_logCPM - (sd_logCPM * sd_penalty))
}

select_reference_set <- function(metrics_tbl, n_ref, prevalence_min, sd_quantile_cut = 0.5) {
  cand <- metrics_tbl %>%
    filter(prevalence >= prevalence_min) %>%
    arrange(desc(median_logCPM))
  
  if (nrow(cand) == 0) stop("No genes passed prevalence filter; lower prevalence_min or check inputs.")
  
  sd_cut <- quantile(cand$sd_logCPM, probs = sd_quantile_cut, na.rm = TRUE)
  cand2 <- cand %>% filter(sd_logCPM <= sd_cut)
  
  if (nrow(cand2) < n_ref) {
    message("Warning: not enough genes after SD filter; relaxing SD filter automatically.")
    cand2 <- cand
  }
  
  cand2 %>%
    arrange(desc(score), desc(median_logCPM)) %>%
    slice_head(n = n_ref)
}

# ----------------------------
# HOST reference set (across ALL samples)
# ----------------------------
host_long <- merged_long %>% filter(genome == host_genome_value)

host_mats <- make_logcpm_matrix(host_long, "GeneID", "sample", "count")
host_metrics <- summarise_gene_stability(host_mats$cpm, host_mats$logcpm,
                                         present_cpm = host_present_cpm,
                                         sd_penalty  = sd_penalty)

host_ref <- select_reference_set(host_metrics, n_ref = n_host,
                                 prevalence_min  = host_prevalence_min,
                                 sd_quantile_cut = sd_quantile_cut)

write_csv(host_metrics, file.path(out_dir, "host_gene_metrics_all_samples.csv"))
write_csv(host_ref,     file.path(out_dir, paste0("host_reference_set_top_", n_host, ".csv")))

# ----------------------------
# PHAGE reference set (INFECTED samples only)
# ----------------------------
phage_long <- merged_long %>%
  filter(genome == phage_genome_value, sample %in% infected_samples)

phage_mats <- make_logcpm_matrix(phage_long, "GeneID", "sample", "count")
phage_metrics <- summarise_gene_stability(phage_mats$cpm, phage_mats$logcpm,
                                          present_cpm = phage_present_cpm,
                                          sd_penalty  = sd_penalty)

phage_ref <- select_reference_set(phage_metrics, n_ref = n_phage,
                                  prevalence_min  = phage_prevalence_min,
                                  sd_quantile_cut = sd_quantile_cut)

write_csv(phage_metrics, file.path(out_dir, "phage_gene_metrics_infected_samples.csv"))
write_csv(phage_ref,     file.path(out_dir, paste0("phage_reference_set_top_", n_phage, ".csv")))

cat("DONE.\n")
cat("Host ref:", file.path(out_dir, paste0('host_reference_set_top_', n_host, '.csv')), "\n")
cat("Phage ref:", file.path(out_dir, paste0('phage_reference_set_top_', n_phage, '.csv')), "\n")
cat("Phage infected samples used:", length(infected_samples), "\n")
