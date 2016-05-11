library(readr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(doMC)
n_chunks = 50
registerDoMC(n_chunks)

line_order = c("A13", "A31", "A41", "A23", "A22", "A42")[6:1]
new_raw_data = read_delim("./data/recalibrated_snps_cohort.table", "\t")

raw_data = new_raw_data %>%
  separate(A13.GT, c("A13_1", "A13_2"), sep = "/") %>%
  separate(A22.GT, c("A22_1", "A22_2"), sep = "/") %>%
  separate(A23.GT, c("A23_1", "A23_2"), sep = "/") %>%
  separate(A31.GT, c("A31_1", "A31_2"), sep = "/") %>%
  separate(A41.GT, c("A41_1", "A41_2"), sep = "/") %>%
  separate(A42.GT, c("A42_1", "A42_2"), sep = "/")

 # raw_data = read_delim("./data/All_biallelic_SNPs.table", " ")
 # raw_data$rows <- NULL

all_poly = raw_data %>%
            select(CHROM, POS, QUAL, TYPE,
                   A13_1, A13_2, A13.GQ,
                   A31_1, A31_2, A31.GQ,
                   A41_1, A41_2, A41.GQ,
                   A23_1, A23_2, A23.GQ,
                   A22_1, A22_2, A22.GQ,
                   A42_1, A42_2, A42.GQ) %>%
            filter(!grepl("_", CHROM)) %>%
            filter(!grepl("chrM", CHROM))

just_snps = all_poly %>% filter(TYPE == "SNP")
not_snps = all_poly %>% filter(TYPE != "SNP")

total = dim(just_snps)[1]
chunk_size = floor(total/n_chunks - 1)
last_chunk = total - (n_chunks - 1) * chunk_size
CHUNKS = data_frame(CHUNKS = factor(c(rep(1:(n_chunks - 1), each = chunk_size), rep(chunk_size, last_chunk))))
just_snps = bind_cols(CHUNKS, just_snps)
n_snp = names(just_snps)
snp_mask = (grepl("_1", n_snp) | grepl("_2", n_snp))
snp_array = as.matrix(just_snps[,snp_mask])
snp_array_list = dlply(just_snps, .(CHUNKS), function(x) as.matrix(x[,snp_mask]))
head(snp_array_list[[1]])

get_biallelic <- function(x){
  tx = table(x)
  if(length(tx) != 2) return(FALSE)
  else if(any(names(tx) == ".")) return(FALSE)
  else if(any(names(tx) == "*")) return(FALSE)
  else return(TRUE)
}

is_bialelic = tbl_df(ldply(snp_array_list,
                    function(snp_array) adply(snp_array, 1, get_biallelic),
                    .parallel = TRUE))

just_snps <- just_snps %>%
     mutate(is_bialelic = is_bialelic$V1) %>%
     select(CHROM, POS, QUAL, TYPE, is_bialelic,
             A13_1, A13_2, A13.GQ,
             A31_1, A31_2, A31.GQ,
             A41_1, A41_2, A41.GQ,
             A23_1, A23_2, A23.GQ,
             A22_1, A22_2, A22.GQ,
             A42_1, A42_2, A42.GQ)

not_snps <- not_snps %>%
     mutate(is_bialelic = FALSE) %>%
     select(CHROM, POS, QUAL, TYPE, is_bialelic,
             A13_1, A13_2, A13.GQ,
             A31_1, A31_2, A31.GQ,
             A41_1, A41_2, A41.GQ,
             A23_1, A23_2, A23.GQ,
             A22_1, A22_2, A22.GQ,
             A42_1, A42_2, A42.GQ)

all_poly = bind_rows(just_snps, not_snps) %>% arrange(CHROM, POS)

bialelic_snps = filter(all_poly, is_bialelic)
total = dim(bialelic_snps)[1]
chunk_size = floor(total/n_chunks - 1)
last_chunk = total - (n_chunks - 1) * chunk_size
CHUNKS = data_frame(CHUNKS = factor(c(rep(1:(n_chunks - 1), each = chunk_size), rep(chunk_size, last_chunk))))
bialelic_snps = bind_cols(CHUNKS, bialelic_snps)
n_snp = names(bialelic_snps)
snp_mask = (grepl("_1", n_snp) | grepl("_2", n_snp))
snp_array = as.matrix(bialelic_snps[,snp_mask])
snp_array_list = dlply(bialelic_snps, .(CHUNKS), function(x) as.matrix(x[,snp_mask]))
head(snp_array_list[[1]])

save(all_poly, bialelic_snps, snp_array_list, file = "./data/just_snps.Rdata")
#load("data/just_snps.Rdata")
