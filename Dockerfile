FROM rocker/r-ver:4.3.2

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    vcftools \
    unzip \
    wget \
    pandoc \
    curl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install PLINK 1.9 manually
RUN curl -sSL https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip -o /tmp/plink.zip && \
    unzip /tmp/plink.zip -d /usr/local/bin && \
    chmod +x /usr/local/bin/plink && \
    rm /tmp/plink.zip

# ✅ Install R packages
RUN R -e "install.packages(c('tidyverse', 'qqman'), repos='https://cloud.r-project.org')"

# Working directory
WORKDIR /home/gwas

# Copy pipeline scripts
COPY scripts/ /home/gwas/scripts/
RUN chmod +x scripts/*.sh

# Create folders
RUN mkdir -p data results

# Run full pipeline by default, then drop into shell
CMD ["bash", "-c", "scripts/gwas_data.sh && scripts/gwas_pipeline.sh && exec bash"]