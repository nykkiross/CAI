# =============================================================================
# Build CAI reference sets from merged_counts CSV
#   - builds a reference set based on ALL counts across all samples for RNAseq
#   - takes consistently highly expressed genes (set to 50, can be adjusted)
# ==============================================================================
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

# if you want an option to also include a phage in analysis (infected vs. uninfected) there is another version of that script that
# includes both host and phage, including extracting infected-only timepoints

# ============
# LOAD COUNTS
# ============

# requires document with all RNAseq counts
merged_counts_csv <- "/path/to/yourfile/merged_counts_ALL.csv"

# if CSV is "wide format" - one column for geneid and one column for each sample - convert to long
wide <- read_csv(merged_counts_csv, show_col_types = FALSE)

# Pivot: all columns except GeneID become sample/count pairs
long <- wide %>%
  pivot_longer(
    cols      = -GeneID,
    names_to  = "sample",
    values_to = "count"
  ) %>%
  mutate(
    GeneID = as.character(GeneID),
    sample = as.character(sample),
    count  = as.numeric(count),
    genome = as.character("your_genome")
    )

gene_col   <- "Geneid"     # e.g., "locus_tag", "gene_id", "feature_id"
sample_col <- "sample"        # e.g., "sample", "well", "sample_id"
count_col  <- "count"         # e.g., "count", "counts", "raw_count"
genome_col <- "genome"        # e.g., "genome", "organism", "species"

organism_genome_value  <- "your_genome"

# ===============
# SET PARAMETERS
# ===============

# Reference set sizes
n_organism  <- 50

# Expression / prevalence thresholds
present_cpm <- 1     # gene considered "present" in a sample if CPM >= this
prevalence_min   <- 0.85  # gene must be present in >=85% samples

# Variability filtering (keep "stable" genes):
# we keep genes with sd_logCPM <= sd_quantile_cut among candidates (lower = stricter)
sd_quantile_cut <- 0.50        # 0.50 = keep lower 50% SD (stable half)

# Scoring: higher median is good; higher SD is bad
sd_penalty <- 0.5              # increase to penalize variability more strongly

# Output folder
out_dir <- "./CAI_reference_sets"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ========
# Helpers
# ========
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

# ==================
# CAI REFERENCE SET
# ==================

mats <- make_logcpm_matrix(long, "GeneID", "sample", "count")
metrics <- summarise_gene_stability(mats$cpm, mats$logcpm,
                                         present_cpm = present_cpm,
                                         sd_penalty  = sd_penalty)

ref_set <- select_reference_set(metrics, n_ref = n_organism,
                                 prevalence_min  = prevalence_min,
                                 sd_quantile_cut = sd_quantile_cut)

write_csv(metrics, file.path(out_dir, "gene_metrics_all_samples.csv"))
write_csv(ref_set,     file.path(out_dir, paste0("reference_set_top_", n_organism, ".csv")))

cat("DONE.\n")
cat("Organism ref:", file.path(out_dir, paste0('reference_set_top_', n_host, '.csv')), "\n")
