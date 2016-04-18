source("./find_private_aleles.R")

by_line = laply(strsplit(colnames(snp_array), "_"), `[`, 1)
get_heterozygote = function(x) laply(tapply(x, by_line, unique), length) - 1

is_heterozygote = aaply(snp_array, 1, get_heterozygote)
colnames(is_heterozygote) <- paste(unique(by_line), "hets", sep = "_")
save(is_heterozygote, file = "./data/heterozygotes_position.Rdata")
load("./data/heterozygotes_position.Rdata")
