# 🧪 CDI Domain Template

Reusable template for all **Complex Data Insights (CDI)** domain repositories, designed for layered Bookdown projects in R and Python.

📘 **Use this template** on GitHub to create new domain-specific guides with consistent structure and automation.

---

## 📦 Features

- ✅ Modular Bookdown structure for:
  - Exploratory Data Analysis (EDA)
  - Visualization (VIZ)
  - Statistical Analysis (STATS)
  - Machine Learning (ML)
- ✅ Ready-to-run R and Python environment setup
- ✅ Preconfigured GitBook and PDF output
- ✅ Clean `.gitignore`, metadata, and styling files
- ✅ Ideal for reproducible data science workflows

---

## 🚀 How to Use This Template

1. Click **"Use this template"** on GitHub.
2. Clone your new repository locally:
   ```bash
   git clone https://github.com/your-org/your-new-domain.git
   cd your-new-domain
   ```
3. Make the setup scripts executable:
   ```bash
   chmod +x scripts/setup_r_env.sh
   chmod +x scripts/setup_py_env.sh
   ```
4. Run the environment setup:
   ```bash
   ./scripts/setup_r_env.sh
   ./scripts/setup_py_env.sh
   ```
5. Build the book:
   ```bash
   ./scripts/build-all.sh
   ```

---

## 🗂️ Folder Overview

```text
data/       → Placeholder dataset (e.g., iris.csv)
images/     → Cover image (use standardized cover.png)
library/    → CSL and BibTeX files
scripts/    → Setup and build automation (R + Python)
*.Rmd       → Modular Q&A content by layer
index-*.Rmd → Bookdown entry files per layer
_output.yml → Bookdown output config
```

---

## 🔧 Notes

- **R packages** managed using `renv` and modular install scripts.
- **Python packages** listed in `requirements.txt`, installed using `venv/`.
- **Cover image** must be named `images/cover.png`.
- Customize `.Rmd` files and `scripts/domain.R` to reflect your specific dataset and goals.

---

## 📄 License

This template is licensed under the [MIT License](LICENSE).  
Attribution appreciated when used for educational or published work.

---

## 🔍 Maintainers

Created and maintained by the **CDI Team** at [ComplexDataInsights.com](https://complexdatainsights.com).
