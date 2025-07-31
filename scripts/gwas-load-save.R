
library(tidyverse)

# Step 1: Load genotype data using readr
ped_data <- read_table("data/sativas413.ped", col_names = FALSE, show_col_types = FALSE)

# Step 2: Save as compressed RDS file
write_rds(ped_data, file = "data/sativas413.rds")

# Step 3: Load metadata files
map_data <- read_table("data/sativas413.map", 
                       col_names = c("chr", "snp_id", "gen_dist", "bp_pos"), 
                       show_col_types = FALSE)

fam_data <- read_table("data/sativas413.fam", 
                       col_names = c("FID", "IID", "PID", "MID", "sex", "phenotype"), 
                       show_col_types = FALSE)

phenotype_data <- read_tsv("data/sativas413_phenotypes.txt", show_col_types = FALSE)