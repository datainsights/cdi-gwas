#!/usr/bin/env Rscript

# 📊 Plot GWAS results using qqman + tidyverse

# --- Load libraries ---
library(tidyverse)
library(qqman)

# --- Parse input arguments ---
args <- commandArgs(trailingOnly = TRUE)
assoc_file <- args[1]
out_dir <- args[2]

# --- Read association results ---
gwas <- read_delim(assoc_file, delim = "\t", show_col_types = FALSE)

# --- Clean and filter (optional: p-value sanity check) ---
gwas_clean <- gwas %>% filter(!is.na(P), P > 0, P < 1)

# --- Manhattan plot ---
manhattan_plot_path <- file.path(out_dir, "manhattan.png")
png(manhattan_plot_path, width = 1000, height = 600)
manhattan(
  gwas_clean,
  main = "GWAS Manhattan Plot",
  col = c("dodgerblue", "gray50"),
  cex = 0.6,
  cex.axis = 0.8
)
dev.off()

# --- QQ plot ---
qq_plot_path <- file.path(out_dir, "qqplot.png")
png(qq_plot_path, width = 800, height = 600)
qq(gwas_clean$P, main = "GWAS QQ Plot")
dev.off()

# --- Done ---
cat("✅ Plots saved to:\n")
cat(" -", manhattan_plot_path, "\n")
cat(" -", qq_plot_path, "\n")