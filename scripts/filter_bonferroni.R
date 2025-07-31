library(tidyverse)

input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

assoc <- read_csv(input_file)
bonf_threshold <- 0.05 / nrow(assoc)
assoc_bonf <- assoc %>% dplyr::filter(P_value < bonf_threshold)

write_csv(assoc_bonf, output_file)