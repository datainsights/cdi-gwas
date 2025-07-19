#!/bin/bash

# Exit immediately if any command fails
set -e

# Print the current R version (useful for logging and debugging)
echo "📦 R version:"
Rscript -e 'R.version.string'

# Function to link the correct Bookdown YAML config file
# Arguments:
#   $1 - The config filename to use (e.g., _bookdown-viz.yml)
link_bookdown_config() {
  local config="$1"
  echo "🔧 Using config: $config"

  # Create or update a symbolic link named _bookdown.yml pointing to the desired config
  cp -f "$config" _bookdown.yml
}

# -------------------------------
# GWAS GitBook
# -------------------------------
if [[ "$1" == "gwas-gitbook" ]]; then
  echo "📘 Building GWAS GitBook..."
  cp -f index-gwas-gitbook.Rmd index.Rmd
  link_bookdown_config _bookdown-gwas.yml
  mkdir -p docs
  Rscript -e 'bookdown::render_book("index.Rmd", "bookdown::gitbook", output_dir = "docs")'
  rm index.Rmd
  echo "✅ GWAS GitBook complete → /docs"


# # -------------------------------
# # GWAS PDF
# # -------------------------------
# elif [[ "$1" == "gwas-pdf" ]]; then
#   echo "📘 Building GWAS PDF..."
#   cp -f index-gwas-pdf.Rmd index.Rmd
#   link_bookdown_config _bookdown-gwas.yml
#   mkdir -p gwas-pdf
#   Rscript -e 'bookdown::render_book("index.Rmd", "bookdown::pdf_book", output_dir = "gwas-pdf")'
#   rm index.Rmd
#   echo "✅ GWAS PDF complete → /gwas-pdf"




# -------------------------------
# Help / fallback
# -------------------------------
else
  echo "❌ Unknown build option: $1"
  echo "Usage: $0 {eda-gitbook|viz-gitbook|viz-pdf|stats-gitbook|stats-pdf|ml-gitbook|ml-pdf}"
  exit 1
fi
