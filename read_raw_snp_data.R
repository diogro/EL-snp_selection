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

just_snps = raw_data %>%
            select(CHROM, POS, A22_hets, A13_1:A42_2) %>%
            filter(!grepl("_", CHROM)) %>%
            filter(!grepl("chrM", CHROM))

snp_array = as.matrix(just_snps[,4:15])

x = snp_array[4,]
is_private = function(x){
  tx = table(x)
  if(any(tx == 2)){
    rare = names(tx)[which(min(tx) == tx)]
    lines = names(x)[x == rare]
    line = unique(laply(strsplit(lines, "_"), `[`, 1))
    if(length(line) == 1) return(c(TRUE, line))
    else return(c(FALSE, NA))
  }
  else return(c(FALSE, NA))
}
#get_private = aaply(snp_array, 1, is_private, .parallel = TRUE)
#save(get_private, file = "./data/private_snp_position.Rdata")
load("./data/private_snp_position.Rdata")

just_snps = mutate(just_snps,
                   is_private = as.logical(get_private[,1]),
                   p_line = get_private[,2])

just_snps %>%
  filter(is_private) %>%
  count(p_line, CHROM) %>%
  ggplot(aes(p_line, log(n), group = p_line)) +
  geom_point() +
  facet_wrap(~CHROM, scales = "free_y")

just_snps %>%
  filter(is_private) %>%
  select(CHROM, POS, p_line) %>%
  filter(CHROM == "chr5") %>%
  ggplot(aes(p_line, POS, color = p_line)) +
  geom_point(size = 0.3) + coord_flip()