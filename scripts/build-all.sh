#!/bin/bash

set -e

echo "🚀 Starting full build for all levels..."

levels=("eda-gitbook" "viz-gitbook" "viz-pdf" "stats-gitbook" "stats-pdf" "ml-gitbook" "ml-pdf")
for level in "${levels[@]}"; do
  echo ""
  echo "🔁 Building: $level"
  bash ./scripts/build-configurable.sh "$level"
done

echo ""
echo "✅ All levels built successfully."

echo "🧹 Cleaning up intermediate build files..."
rm -f inter-pdf/*.md inter-pdf/*.tex
rm -f adv-pdf/*.md adv-pdf/*.tex
rm -rf inter-pdf/_bookdown_files/
rm -rf adv-pdf/_bookdown_files/
rm -f ./*.rds
echo "✅ Cleanup complete."

# Remove symlink to avoid accidental reuse
rm -f index.Rmd
rm -f _bookdown.yml