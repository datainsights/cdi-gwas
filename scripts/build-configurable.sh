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
# EDA GitBook
# -------------------------------
if [[ "$1" == "eda-gitbook" ]]; then
  echo "📘 Building EDA GitBook..."
  cp -f index-eda-gitbook.Rmd index.Rmd
  link_bookdown_config _bookdown-eda.yml
  mkdir -p docs
  Rscript -e 'bookdown::render_book("index.Rmd", "bookdown::gitbook", output_dir = "docs")'
  rm index.Rmd
  echo "✅ EDA GitBook complete → /docs"

# -------------------------------
# VIZ GitBook
# -------------------------------
elif [[ "$1" == "viz-gitbook" ]]; then
  echo "📘 Building Visualization GitBook..."
  cp -f index-viz-gitbook.Rmd index.Rmd
  link_bookdown_config _bookdown-viz.yml
  mkdir -p viz-gitbook
  Rscript -e 'bookdown::render_book("index.Rmd", "bookdown::gitbook", output_dir = "viz-gitbook")'
  rm index.Rmd
  echo "✅ Visualization GitBook complete → /viz-gitbook"

# -------------------------------
# Visualization PDF
# -------------------------------
elif [[ "$1" == "viz-pdf" ]]; then
  echo "📘 Building Visualization PDF..."
  cp -f index-viz-pdf.Rmd index.Rmd
  link_bookdown_config _bookdown-viz.yml
  mkdir -p viz-pdf
  Rscript -e 'bookdown::render_book("index.Rmd", "bookdown::pdf_book", output_dir = "viz-pdf")'
  rm index.Rmd
  echo "✅ Visualization PDF complete → /viz-pdf"

# -------------------------------
# Statistical Analysis GitBook
# -------------------------------
elif [[ "$1" == "stats-gitbook" ]]; then
  echo "📘 Building Statistical Analysis GitBook..."
  cp -f index-stats-gitbook.Rmd index.Rmd
  link_bookdown_config _bookdown-stats.yml
  mkdir -p stats-gitbook
  Rscript -e 'bookdown::render_book("index.Rmd", "bookdown::gitbook", output_dir = "stats-gitbook")'
  rm index.Rmd
  echo "✅ Statistical Analysis GitBook complete → /stats-gitbook"

# -------------------------------
# Statistical Analysis PDF
# -------------------------------
elif [[ "$1" == "stats-pdf" ]]; then
  echo "📘 Building Statistical Analysis PDF..."
  cp -f index-stats-pdf.Rmd index.Rmd
  link_bookdown_config _bookdown-stats.yml
  mkdir -p stats-pdf
  Rscript -e 'bookdown::render_book("index.Rmd", "bookdown::pdf_book", output_dir = "stats-pdf")'
  rm index.Rmd
  echo "✅ Statistical Analysis PDF complete → /stats-pdf"


# -------------------------------
# Machine Learning GitBook
# -------------------------------
elif [[ "$1" == "ml-gitbook" ]]; then
  echo "📘 Building Machine Learning GitBook..."
  cp -f index-ml-gitbook.Rmd index.Rmd
  link_bookdown_config _bookdown-ml.yml
  mkdir -p ml-gitbook
  Rscript -e 'bookdown::render_book("index.Rmd", "bookdown::gitbook", output_dir = "ml-gitbook")'
  rm index.Rmd
  echo "✅ Machine Learning GitBook complete → /ml-gitbook"

# -------------------------------
# Machine Learning PDF
# -------------------------------
elif [[ "$1" == "ml-pdf" ]]; then
  echo "📘 Building Machine Learning PDF..."
  cp -f index-ml-pdf.Rmd index.Rmd
  link_bookdown_config _bookdown-ml.yml
  mkdir -p ml-pdf
  Rscript -e 'bookdown::render_book("index.Rmd", "bookdown::pdf_book", output_dir = "ml-pdf")'
  rm index.Rmd
  echo "✅ Machine Learning PDF complete → /ml-pdf"




# -------------------------------
# Help / fallback
# -------------------------------
else
  echo "❌ Unknown build option: $1"
  echo "Usage: $0 {eda-gitbook|viz-gitbook|viz-pdf|stats-gitbook|stats-pdf|ml-gitbook|ml-pdf}"
  exit 1
fi

echo "🧹 Cleaning up Visualization build files..."
rm -f viz-pdf/*.md viz-pdf/*.tex
rm -f stats-pdf/*.md stats-pdf/*.tex
rm -f ml-pdf/*.md ml-pdf/*.tex
rm -rf viz-pdf/_bookdown_files/
rm -rf stats-pdf/_bookdown_files/
rm -rf ml-pdf/_bookdown_files/
rm -f ./*.rds
echo "✅ Cleanup complete."

# Remove symlink to avoid accidental reuse
rm -f index.Rmd
rm -f _bookdown.yml