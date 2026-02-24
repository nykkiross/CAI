# CAI
generate codon tables and reference gene sets, calculate CAI - ONE organism

This repository contains 3 scripts - one to generate expression-weighted codon tables for ONE organism (using CPM from raw counts data), one to generate a reference set of genes from RNAseq expression data, and one to generate a codon table from the reference genes and calculate Codon Adaptation Index (CAI) for expression-weighted codon usage

Can also be adapted for host/phage pairs - CAI_refsets.R and CAI_calc.R assume phage and host pairs, but phage can be removed if only one organism is being investigated

Determines how well-adapted the actual codons used are to the codon usage of the most highly-expressed genes in the genome

Uses expression data for accuracy in expression

If you have any questions or issues, feel free to contact me at nicole.ross@ufl.edu :) 
