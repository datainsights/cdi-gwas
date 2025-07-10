# 📘 Domain-Specific R Packages for Genome-Wide Association Studies (GWAS)

# Load or install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
library(BiocManager)

# Load or install renv
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
library(renv)

# --- Bioconductor packages for GWAS ---
bioc_pkgs <- c(
  "SNPRelate", 
  "gdsfmt"
  )

# --- CRAN packages for GWAS ---
cran_pkgs <- c(
  "rrBLUP", 
  "BGLR", "DT", 
  "dplyr", 
  "qqman", 
  "poolr", 
  "glue"
  )

# Install Bioconductor packages
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(glue::glue("Installing Bioconductor package: {pkg}"))
    renv::install(pkg, repos = BiocManager::repositories())
  } else {
    message(glue::glue("✔ {pkg} already installed."))
  }
}

# Install CRAN packages
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(glue::glue("Installing CRAN package: {pkg}"))
    renv::install(pkg)
  } else {
    message(glue::glue("✔ {pkg} already installed."))
  }
}

# --- Optional GitHub packages ---
# if (!requireNamespace("pkgname", quietly = TRUE)) {
#   message("Installing GitHub package: username/pkgname")
#   renv::install("username/pkgname")
# }

message("✅ GWAS domain package setup complete.")