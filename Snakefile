# Snakefile

configfile: "config/config.yaml"
prefix = config["prefix"]

rule all:
    input:
        "results/significant_snps_bonferroni.csv",
        "results/significant_snps.png"

rule filter_bonferroni:
    input:
        assoc="data/gwas_results.csv"
    output:
        filtered="results/significant_snps_bonferroni.csv"
    script:
        "scripts/filter_bonferroni.R"

# rule plot_snps:
#     input:
#         csv="results/significant_snps_bonferroni.csv"
#     output:
#         fig="results/significant_snps.png"
#     script:
#         "scripts/plot_snps.R"