# Load necessary library
library(tidyverse)

# Step 1: Load genotype and SNP info (if not already in memory)
ped_data <- read_rds("data/sativas413.rds")
snp_info <- read_table("data/sativas413.map", 
                       col_names = c("chr", "snp_id", "gen_dist", "bp_pos"), 
                       show_col_types = FALSE)

# Step 2: Separate sample IDs and genotype calls
sample_ids <- ped_data[, 1:2]              # FID and IID
genotype_matrix <- ped_data[, -(1:6)]      # Alleles only

# Step 3: Verify expected SNP count
n_snps <- ncol(genotype_matrix) / 2
stopifnot(n_snps == nrow(snp_info))

# Step 4: Combine each pair of allele columns into genotype strings
genotype_calls <- map_dfc(seq(1, ncol(genotype_matrix), by = 2), function(i) {
  paste(genotype_matrix[[i]], genotype_matrix[[i + 1]])
})
names(genotype_calls) <- snp_info$snp_id  # Use SNP IDs as column names

# Step 5: Combine with sample IDs
genotype_tidy <- bind_cols(sample_ids, genotype_calls)


# Step 6: Drop FID and IID from genotype_tidy to isolate genotype columns
geno_alleles <- genotype_tidy[, -c(1, 2)]

# Step 7: Convert allele strings to numeric minor allele counts
geno_minor_allele_count <- map_dfc(geno_alleles, function(allele_vec) {
  # Split all genotype strings (e.g., "A G") into individual alleles
  alleles <- unlist(str_split(allele_vec, " "))
  allele_counts <- table(alleles)

  # Skip SNPs that are monomorphic or malformed
  if (length(allele_counts) < 2) return(rep(NA, length(allele_vec)))

  # Identify the minor allele (less frequent)
  minor_allele <- names(sort(allele_counts))[1]

  # Count how many copies of the minor allele are in each genotype
  sapply(allele_vec, function(gt) {
    if (gt %in% c("0 0", "0 1", "1 0", "1 1", "0", "1")) return(NA)  # filter malformed
    split_alleles <- unlist(str_split(gt, " "))
    if (length(split_alleles) != 2) return(NA)
    sum(split_alleles == minor_allele)
  })
})

# Step 8: Add back sample identifiers
genotype_count <- bind_cols(genotype_tidy[, 1:2], geno_minor_allele_count)

# Step 9: Remove sample columns (FID, IID)
count_only <- genotype_count[, -c(1, 2)]

# Step 10: Filter SNPs by missingness (e.g., keep SNPs with <10% missing values)
snp_missing <- colMeans(is.na(count_only))
snp_keep <- names(snp_missing[snp_missing < 0.1])
filtered_count <- count_only[, snp_keep]

# Step 11: Filter SNPs by minor allele frequency (MAF >= 0.05)
calc_maf <- function(x) {
  p <- mean(x, na.rm = TRUE) / 2
  min(p, 1 - p)
}
snp_maf <- map_dbl(filtered_count, calc_maf)
maf_keep <- names(snp_maf[snp_maf >= 0.05])
final_count <- filtered_count[, maf_keep]

# Step 12: Reattach FID and IID
filtered_geno <- bind_cols(genotype_count[, 1:2], final_count)

# Step 13: Summary of filtering
cat("Original SNPs:", ncol(count_only), "\n")
cat("After missing filter:", length(snp_keep), "\n")
cat("After MAF filter:", length(maf_keep), "\n")

# Step 14: Extract dosage matrix (without FID/IID)
dosage_matrix <- filtered_geno[, -c(1, 2)]

# Step 15: Impute missing values using column means
imputed_matrix <- dosage_matrix %>%
  mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

# Step 16: Add back FID and IID
geno_imputed <- bind_cols(filtered_geno[, 1:2], imputed_matrix)


# Step 17: Extract genotype matrix (exclude FID and IID)
geno_numeric <- geno_imputed[, -c(1, 2)]

# Step 18: Perform PCA using prcomp
pca_result <- prcomp(geno_numeric, center = TRUE, scale. = TRUE)

# Step 19: Combine first 5 PCs with sample IDs
pca_df <- geno_imputed[, 1:2] %>%  # FID and IID
  bind_cols(as_tibble(pca_result$x[, 1:5]))  # PC1 to PC5

# Step 20: Plot PC1 vs PC2
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "PCA of Genotype Data", x = "PC1", y = "PC2") +
  theme_minimal()

# Step 21: Load and align phenotype data with sample metadata (by FID)
fam_data <- read_table("data/sativas413.fam", 
                       col_names = c("FID", "IID", "PID", "MID", "sex", "phenotype"), 
                       show_col_types = FALSE)

phenotype_data <- read_tsv("data/sativas413_phenotypes.txt", show_col_types = FALSE) %>%
  rename(FID = 1)  # Sample IDs in phenotype file match fam_data$FID

sample_metadata <- left_join(fam_data, phenotype_data, by = "FID")

# Step 22: Rename PCA columns to include proper IDs
pca_df <- pca_df %>%
  rename(FID = 1, IID = 2)

# Step 23: Merge PCA data into sample metadata
sample_data <- left_join(sample_metadata, pca_df, by = c("FID", "IID"))

# Step 24: Standardize ID columns in genotype data
geno_imputed <- geno_imputed %>%
  rename(FID = 1, IID = 2)

# Step 25: Merge genotype with metadata
geno_data <- geno_imputed[, -1]  # Drop FID, keep IID and SNPs
gwas_input <- left_join(sample_data, geno_data, by = "IID")

# Step 26: Select trait and covariates
trait <- "Plant height"  # Use the column name as a string
covariates <- c("PC1", "PC2", "PC3")

# Save the merged GWAS input to an RDS file for future use
saveRDS(gwas_input, file = "data/gwas_input.rds")

# Step 27: Construct the model formula
snp_name <- names(geno_data)[2]  # Replace with desired SNP
formula_str <- paste0("`", trait, "` ~ ", snp_name, " + ", paste(covariates, collapse = " + "))
model <- lm(as.formula(formula_str), data = gwas_input)

# Step 28: View model summary
summary(model)


# Perform a GWAS scan
# Step 1: Identify SNP columns using prefix pattern
snp_cols <- grep("^(id|ud|wd|dd|fd)[0-9]+$", names(gwas_input), value = TRUE)

# Step 2: Check number of SNPs selected
length(snp_cols)      # Should return 3755
head(snp_cols, 5)     # Preview first 5 SNPs

# Step 3: Initialize list to collect GWAS results
gwas_results <- list()

# Step 4: Loop over each SNP and fit linear model
for (snp in snp_cols) {
  # Construct formula dynamically
  formula <- as.formula(paste("`Plant height` ~ PC1 + PC2 + PC3 +", snp))
  
  # Fit model safely
  model <- tryCatch(
    lm(formula, data = gwas_input),
    error = function(e) NULL
  )
  
  # Step 5: If successful, extract coefficient and p-value
  if (!is.null(model)) {
    coef_table <- summary(model)$coefficients
    snp_row <- tail(rownames(coef_table), 1)
    
    gwas_results[[snp]] <- tibble(
      SNP = snp,
      Estimate = coef_table[snp_row, "Estimate"],
      P_value = coef_table[snp_row, "Pr(>|t|)"]
    )
  }
}

# Step 6: Combine results and sort by p-value
gwas_df <- bind_rows(gwas_results) %>%
  arrange(P_value)

# Step 7: Save the GWAS results for visualization
write_csv(gwas_df, "data/gwas_results.csv")

# Manhattan Plot using ggplot2
# Step 1: Load GWAS results and SNP position data
gwas_df <- read_csv("data/gwas_results.csv")
map_df <- read_tsv("data/sativas413.map", 
                   col_names = c("CHR", "SNP", "GEN_DIST", "BP_POS"),
                   show_col_types = FALSE)

# Step 2: Merge GWAS results with chromosome position
gwas_annotated <- left_join(gwas_df, map_df, by = "SNP") %>%
  drop_na()  # Remove SNPs with missing position

# Step 3: Compute cumulative position for plotting across chromosomes
gwas_annotated <- gwas_annotated %>%
  arrange(CHR, BP_POS) %>%
  group_by(CHR) %>%
  mutate(BP_CUM = BP_POS + ifelse(row_number() == 1, 0, lag(cumsum(BP_POS), default = 0))) %>%
  ungroup()

# Step 4: Compute –log10(p-value) and color group
gwas_annotated <- gwas_annotated %>%
  mutate(logP = -log10(P_value),
         CHR = as.factor(CHR),
         color_group = as.integer(CHR) %% 2)

# Step 5: Plot Manhattan plot
ggplot(gwas_annotated, aes(x = BP_CUM, y = logP, color = as.factor(color_group))) +
  geom_point(alpha = 0.7, size = 1.2) +
  scale_color_manual(values = c("#003b4a", "dodgerblue")) +
  labs(title = "Manhattan Plot Using ggplot2",
       x = "Genomic Position", y = expression(-log[10](p))) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

# Manhattan Plot using qqman
# Step 1: Load GWAS results
gwas_df <- read_csv("data/gwas_results.csv")

# Step 2: Load SNP position data (.map file)
map_df <- read_tsv("data/sativas413.map", 
                   col_names = c("CHR", "SNP", "GEN_DIST", "BP"),
                   show_col_types = FALSE)

# Step 3: Merge results with map data
gwas_annotated <- left_join(gwas_df, map_df, by = "SNP") %>%
  select(SNP, CHR, BP, P = P_value) %>%
  drop_na()

# Step 4: Create Manhattan plot using qqman
manhattan(gwas_annotated,
          main = "Manhattan Plot of GWAS Results",
          col = c("grey30", "dodgerblue"),
          cex = 0.6,
          cex.axis = 0.9,
          las = 1)

# QQ Plot by qqman
# Load required libraries
library(tidyverse)
library(qqman)

# Step 1: Load GWAS results
gwas_df <- read_csv("data/gwas_results.csv")

# Step 2: Create QQ plot using qqman
qq(gwas_df$P_value,
   main = "QQ Plot of GWAS Results (qqman)")

# QQ plot by ggplo2
# Step 1: Load GWAS results
gwas_df <- read_csv("data/gwas_results.csv")

# Step 2: Calculate expected vs observed -log10(p)
gwas_df <- gwas_df %>%
  filter(!is.na(P_value)) %>%
  mutate(observed = -log10(sort(P_value)),
         expected = -log10(ppoints(n())))

# Step 3: Create QQ plot with reference line using ggplot2
ggplot(gwas_df, aes(x = expected, y = observed)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_point(size = 1.2, alpha = 0.6, color = "steelblue") +
  labs(title = "QQ Plot of GWAS Results (ggplot2)",
       x = "Expected -log10(p)",
       y = "Observed -log10(p)") +
  theme_minimal(base_size = 14)

  # Multiple GWAS correction tests
  # Load required packages
library(tidyverse)

# Step 1: Load GWAS results
gwas_df <- read_csv("data/gwas_results.csv")

# Step 2: Add Bonferroni-corrected threshold
n_tests <- nrow(gwas_df)
bonf_threshold <- 0.05 / n_tests

# Step 3: Apply FDR correction using p.adjust
gwas_df <- gwas_df %>%
  mutate(FDR = p.adjust(P_value, method = "BH"))

# Step 4: Extract significant SNPs
significant_bonf <- gwas_df %>%
  filter(P_value < bonf_threshold)

significant_fdr <- gwas_df %>%
  filter(FDR < 0.05)

# Step 5: Output summary
cat("Bonferroni threshold:", bonf_threshold, "\n")
cat("Number of SNPs passing Bonferroni:", nrow(significant_bonf), "\n")
cat("Number of SNPs passing FDR < 0.05:", nrow(significant_fdr), "\n")

# Volacano plot by ggplot2
# Step 1: Load GWAS results
gwas_df <- read_csv("data/gwas_results.csv")

# Step 2: Compute –log10(p-value)
gwas_df <- gwas_df %>%
  mutate(logP = -log10(P_value))

# Step 3: Create volcano plot
ggplot(gwas_df, aes(x = Estimate, y = logP)) +
  geom_point(alpha = 0.6, color = "grey40") +
  geom_hline(yintercept = -log10(0.05 / nrow(gwas_df)), linetype = "dashed", color = "red") +  # Bonferroni line
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +  # Null effect line
  labs(title = "Volcano Plot of GWAS Results",
       x = "Effect Size (Estimate)",
       y = expression(-log[10](p))) +
  theme_minimal(base_size = 14)


# Add significance status
gwas_df <- gwas_df %>%
  mutate(significant = P_value < 0.05 / nrow(gwas_df))

# Re-plot with color by significance
ggplot(gwas_df, aes(x = Estimate, y = logP, color = significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("grey70", "red")) +
  geom_hline(yintercept = -log10(0.05 / nrow(gwas_df)), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  labs(title = "Volcano Plot with Bonferroni Threshold",
       x = "Effect Size (Estimate)",
       y = expression(-log[10](p))) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank())

  # Significant SNP hits
  # Step 1: Load GWAS results
gwas_df <- read_csv("data/gwas_results.csv")

# Step 2: Calculate Bonferroni threshold
n_tests <- nrow(gwas_df)
bonf_threshold <- 0.05 / n_tests  # ≈ 1.33e-05 for 3755 SNPs

# Step 3: Filter significant hits
significant_hits <- gwas_df %>%
  filter(P_value < bonf_threshold) %>%
  arrange(P_value)

# Step 4: Load SNP position data
map_df <- read_tsv("data/sativas413.map",
                   col_names = c("CHR", "SNP", "GEN_DIST", "BP"),
                   show_col_types = FALSE)

# Step 5: Annotate significant SNPs with position
annotated_hits <- left_join(significant_hits, map_df, by = "SNP") %>%
  select(SNP, CHR, BP, Estimate, P_value)

# Step 6: Save for downstream analysis
write_csv(annotated_hits, "data/significant_snps_bonferroni.csv")
