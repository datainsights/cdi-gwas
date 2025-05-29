# 🧬 CDI: Genome-Wide Association Studies (GWAS)

[![Live Site](https://img.shields.io/badge/visit-site-blue?logo=githubpages)](https://datainsights.github.io/cdi-gwas)

📘 **Live guide:** [https://datainsights.github.io/cdi-gwas](https://datainsights.github.io/cdi-gwas)

Modular, reproducible framework for GWAS workflows using  
**Bookdown**, **renv**, and **GitHub Actions** — deployed via GitHub Pages.

---

![Test Book Build](https://github.com/datainsights/cdi-gwas/actions/workflows/test-book.yml/badge.svg)
![Deploy Book](https://github.com/datainsights/cdi-gwas/actions/workflows/deploy-book.yml/badge.svg)

> Manual test + auto-deploy CI split for clean and controlled builds.

---

## 📘 Overview

This domain guides users through Genome-Wide Association Studies (GWAS) using structured layers and reproducible analysis.

Layered structure:

- 🔍 **Exploratory Data Analysis (EDA)** layer
- 📊 **Visualization (VIZ)** layer *(coming soon)*
- 📐 **Statistical Analysis (STATS)** layer *(coming soon)*
- 🧠 **Machine Learning (ML)** layer *(coming soon)*

---

## 🛠️ Environment Setup

This project supports both **R** and **Python** workflows.

### 🔄 Option 1: Restore R dependencies directly

••Requires existing renv.lock••

```bash
Rscript -e 'renv::restore()'
```
### ⚙️ Option 2: Run the full environment setup (recommended)

```bash
# Make setup scripts executable
chmod +x scripts/setup_r_env.sh
chmod +x scripts/setup_py_env.sh

# Run R and Python environment setup
./scripts/setup_r_env.sh
./scripts/setup_py_env.sh
```

### 📦 Notes

- **R packages** are managed with `renv` and modular installer scripts in `scripts/`.
- **Python packages** are listed in `requirements.txt` and installed using a virtual environment (`venv/`) created with Python’s built-in `venv` module.
- Customize:
  - `scripts/common.R` and `scripts/gwas.R` (for R)
  - `requirements.txt` (for Python)

---

## 📁 Data Sources

Datasets will include typical GWAS inputs: genotype matrices, phenotype data, and annotations. Place them in the `data/` folder.

---

## 📄 License

This guide is released under the [MIT License](LICENSE).  
Attribution appreciated when used in educational or published works.

---

## 🧠 Maintainers

Created and maintained by the **CDI Team** at [ComplexDataInsights.com](https://complexdatainsights.com)
