args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

library(tidyverse)
library(ggrepel)

snps <- read_csv(input_file)
snps <- snps %>% mutate(logP = -log10(P_value))

ggplot(snps, aes(x = SNP, y = logP)) +
  geom_point(color = "#2C7BB6", size = 3) +
  geom_text_repel(data = snps %>% filter(logP > 12), aes(label = SNP), size = 3) +
  labs(title = "Significant SNPs (Bonferroni)",
       x = "SNP",
       y = expression(-log[10](P))) +
  theme_minimal()

ggsave(output_file, width = 8, height = 5)