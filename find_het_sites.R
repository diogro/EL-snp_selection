source("./find_private_aleles.R")

by_line = laply(strsplit(colnames(snp_array), "_"), `[`, 1)
get_heterozygote = function(x) tapply(x, by_line, function(row) length(unique(row))) - 1

is_heterozygote = tbl_df(ldply(snp_array_list,
                               function(snp_array) aaply(snp_array, 1, get_heterozygote),
                               .parallel = TRUE))
names(is_heterozygote) <- paste(unique(by_line), "hets", sep = "_")
save(is_heterozygote, file = "./data/heterozygotes_position.Rdata")
#load("./data/heterozygotes_position.Rdata")

just_snps <- bind_cols(just_snps, is_heterozygote)

just_hets <- just_snps %>% select(CHROM, POS, contains("_hets"))
just_hets

library(RcppRoll)
current_chr = "chr1"
mean_het <-
  just_hets %>%
  filter(CHROM == current_chr) %>%
  select(contains("_hets")) %>%
  {roll_mean(as.matrix(.), 100001,
             fill = NA, align = "center")}
names(mean_het) <- paste(unique(by_line), "mean_het", sep = "_")
