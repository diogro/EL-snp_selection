library(readr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(doMC)
registerDoMC(21)

# new_raw_data = read_delim("./data/recalibrated_t90_snps_cohort.table", "\t")
#
# raw_data = new_raw_data %>%
#   separate(A13.GT, c("A13_1", "A13_2"), sep = "/") %>%
#   separate(A22.GT, c("A22_1", "A22_2"), sep = "/") %>%
#   separate(A23.GT, c("A23_1", "A23_2"), sep = "/") %>%
#   separate(A31.GT, c("A31_1", "A31_2"), sep = "/") %>%
#   separate(A41.GT, c("A41_1", "A41_2"), sep = "/") %>%
#   separate(A42.GT, c("A42_1", "A42_2"), sep = "/")
#
# # raw_data = read_delim("./data/All_biallelic_SNPs.table", " ")
# # raw_data$rows <- NULL
#
line_order = c("A13", "A31", "A41", "A23", "A22", "A42")[6:1]
#
# just_snps = raw_data %>%
#             select(CHROM, POS, A13_1:A42_2) %>%
#             filter(!grepl("_", CHROM)) %>%
#             filter(!grepl("chrM", CHROM))
#
# snp_array = as.matrix(just_snps[,3:14])
# snp_array_list = dlply(just_snps, .(CHROM), function(x) as.matrix(x[,3:14]))
#
# x = snp_array_list[[19]][1,]
# get_biallelic <- function(x){
#   tx = table(x)
#   if(length(tx) != 2) return(FALSE)
#   else if(any(names(tx) == ".")) return(FALSE)
#   else return(TRUE)
# }
#
# is_bialelic = ldply(snp_array_list,
#                     function(snp_array) adply(snp_array, 1, get_biallelic),
#                     .parallel = TRUE)
# bind_cols(just_snps, select(is_bialelic, V1))
#
# just_snps = bind_cols(just_snps, select(is_bialelic, V1)) %>%
#   filter(V1) %>%
#   select(CHROM, POS, A13_1:A42_2)
#
# snp_array = as.matrix(just_snps[,3:14])
# snp_array_list = dlply(just_snps, .(CHROM), function(x) as.matrix(x[,3:14]))
#
# save(just_snps, snp_array_list, file = "./data/just_snps.Rdata")
load("data/just_snps.Rdata")
