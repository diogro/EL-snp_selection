library(readr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(doMC)
registerDoMC(4)

raw_data = read_delim("./data/All_biallelic_SNPs.table", " ")
raw_data$rows <- NULL

line_order = c("A13", "A31", "A41", "A23", "A22", "A42")[6:1]

just_snps = raw_data %>%
            select(CHROM, POS, A13_1:A42_2) %>%
            filter(!grepl("_", CHROM)) %>%
            filter(!grepl("chrM", CHROM))

snp_array = as.matrix(just_snps[,3:14])
snp_array_list = dlply(just_snps, .(CHROM), function(x) as.matrix(x[,3:14]))