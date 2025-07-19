#!/bin/bash

# 🚀 Prepare rice GWAS genotype and phenotype data (PLINK format)

# --- Paths and filenames ---
ZIP_URL="http://ricediversity.org/data/sets/44kgwas/RiceDiversity.44K.MSU6.Genotypes_PLINK.zip"
ZIP_FILE="data/rice_gwas_genotypes.zip"
EXTRACT_DIR="data/RiceDiversity_44K_Genotypes_PLINK"
PHENO_URL="http://www.ricediversity.org/data/sets/44kgwas/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt"
PHENO_OUT="data/sativa413_phenotypes.txt"

# --- Step 1: Ensure data folder exists ---
mkdir -p data

# --- Step 2: Download genotype zip file (if not already present) ---
if [ ! -f "$ZIP_FILE" ]; then
  echo "⬇️ Downloading genotype data..."
  wget --no-check-certificate -O "$ZIP_FILE" "$ZIP_URL"
else
  echo "✅ Genotype zip already exists: $ZIP_FILE"
fi

# --- Step 3: Unzip genotype data ---
echo "📂 Extracting genotype files..."
unzip -o "$ZIP_FILE" -d data/

# --- Step 4: Move files up and clean nested folder ---
if [ -d "$EXTRACT_DIR" ]; then
  mv "$EXTRACT_DIR"/* data/
  rm -rf "$EXTRACT_DIR"
fi
rm -rf data/__MACOSX
rm -f "$ZIP_FILE"

# --- Step 5: Download phenotype file and rename ---
echo "⬇️ Downloading phenotype file..."
wget --no-check-certificate -P data/ "$PHENO_URL"
mv data/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt "$PHENO_OUT"

echo "✅ GWAS data successfully prepared in the data/ folder."