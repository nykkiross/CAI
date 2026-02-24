# ================
# CAI CALCULATION
# ================

# host using host ref set
# phage using host ref set
# phage using phage ref set

# ==================
# CODON TO AA TABLE
# ==================
library(Biostrings)
library(tibble)

codon_table <- tibble(
  codon = tolower(names(Biostrings::GENETIC_CODE)),
  aa    = Biostrings::GENETIC_CODE
)
# filter out stop codons
codon_table <- codon_table %>%
  filter(aa != "*")

codon_counts_weighted <- codon_counts_weighted %>%
  left_join(codon_table, by = "codon")

# ==================================
# GENERATE CODON TABLES FOR REF SET
# ==================================
host_ref_ids <- sub("^gene-", "", host_ref$gene)

phage_ref_ids <- phage_ref %>%
  transmute(
    phage_gp = gene %>%
      sub("^gene-", "", .) %>%         # remove leading "gene-"
      sub("^phage_locus_prefix_", "", .)            # remove leading "prefix_" if necessary; can be skipped
  )

host_ref_codon <- codon_counts_weighted %>%
  filter(genome == "host", locus_tag_join %in% host_ref_ids) %>%
  mutate(weighted_count = count * CPM) %>%
  group_by(codon) %>%
  summarise(count = sum(weighted_count, na.rm = TRUE), .groups = "drop") %>%
  left_join(codon_table, by = "codon") %>%
  filter(!is.na(aa), aa !="*") %>%
  group_by(codon, aa) %>%
  summarise(count = sum(count), .groups = "drop")

phage_ref_codon <- codon_counts_weighted %>%
  filter(genome == "phage") %>%
  mutate(
    phage_gp = locus_tag_join %>% sub("^gene-", "", .)  # just in case
  ) %>%
  semi_join(phage_ref_ids, by = "phage_gp") %>%
  mutate(weighted_count = count * CPM) %>%
  group_by(codon) %>%
  summarise(count = sum(weighted_count, na.rm = TRUE), .groups = "drop") %>%
  left_join(codon_table, by = "codon") %>%
  filter(!is.na(aa), aa != "*")

# ===============================================
# CONVERT TO RELATIVE ADAPTIVENESS WEIGHTS (w_i)
# ===============================================
host_weights <- host_ref_codon %>%
  group_by(aa, .drop = TRUE) %>%     # <— key
  mutate(
    freq = count / sum(count),
    w    = freq / max(freq)
  ) %>%
  ungroup()

phage_weights <- phage_ref_codon %>%
  group_by(aa) %>%
  mutate(
    freq = count / sum(count),
    w    = freq / max(freq)
  ) %>%
  ungroup()

# =========================================================
# CALCULATE CAI PER GENE USING codon_counts_weighted TABLE
# =========================================================
# function to compute
calc_cai <- function(codon_df, weights_df) {
  
  codon_df %>%
    left_join(weights_df %>% select(codon, w), by = "codon") %>%
    mutate(w = ifelse(is.na(w), 0.01, w)) %>%  # avoid log(0)
    group_by(genome, locus_tag_join) %>%
    summarise(
      CAI = exp(mean(log(w))),
      .groups = "drop"
    )
}

# run on host using host ref
host_cai <- calc_cai(
  codon_counts_weighted %>% filter(genome == "2457T"),
  host_weights
)

# run on phage using host ref
phage_host_cai <- calc_cai(
  codon_counts_weighted %>% filter(genome == "Sf14"),
  host_weights
)

# run on phage using phage ref
phage_phage_cai <- calc_cai(
  codon_counts_weighted %>% filter(genome == "Sf14"),
  phage_weights
)

write.csv(host_cai, "/Users/nross1994/Documents/Doore Lab/RNAseq Analysis/host_CAI.csv")
write.csv(phage_host_cai, "/Users/nross1994/Documents/Doore Lab/RNAseq Analysis/phage_CAI_hostref.csv")
write.csv(phage_phage_cai, "/Users/nross1994/Documents/Doore Lab/RNAseq Analysis/phage_CAI_phageref.csv")
