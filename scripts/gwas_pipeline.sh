#!/bin/bash

# GWAS Pipeline - CLI Version
# Author: TMB (CDI)
# Description: Run a basic GWAS analysis using PLINK and generate summary plots in R.

set -e
set -o pipefail

echo "📦 Starting GWAS pipeline..."
echo "🕒 Timestamp: $(date)"
echo

# ------------------------------------------
# 1. Define Input Paths
# ------------------------------------------
DATA_DIR="data"
OUT_DIR="results"
mkdir -p "$OUT_DIR"

# Input file paths (PED/MAP format)
GENOTYPE="$DATA_DIR/sativas413"                 # Expects .ped/.map/.fam
PHENOTYPE="$DATA_DIR/sativa413_phenotypes.txt" # Tab-delimited with FID, IID, PHENO
COVARIATE="$DATA_DIR/covariates.txt"           # Optional

# ------------------------------------------
# 2. Run GWAS using PLINK (direct PED/MAP input)
# ------------------------------------------
echo "🔍 Running association test with PLINK (PED/MAP format)..."

plink \
  --file "$GENOTYPE" \
  --pheno "$PHENOTYPE" \
  --pheno-name PHENO \
  --covar "$COVARIATE" \
  --covar-name Age,Sex,PC1,PC2,PC3 \
  --logistic \
  --allow-no-sex \
  --out "$OUT_DIR/gwas_results"

echo "✅ PLINK GWAS complete."

# ------------------------------------------
# 3. Generate Manhattan & QQ Plots in R
# ------------------------------------------
echo "📊 Generating Manhattan and QQ plots..."

Rscript scripts/gwas_plots.R "$OUT_DIR/gwas_results.assoc.logistic" "$OUT_DIR"

echo "✅ Plots saved to $OUT_DIR"
echo "🎉 GWAS pipeline completed!"