source("./find_private_aleles.R")

by_line = laply(strsplit(colnames(snp_array), "_"), `[`, 1)
get_heterozygote = function(x) tapply(x, by_line, function(row) length(unique(row))) - 1

is_heterozygote = ldply(snp_array_list,
                        function(snp_array) aaply(snp_array, 1, get_heterozygote),
                        .parallel = TRUE)[,-1]
names(is_heterozygote) <- paste(unique(by_line), "hets", sep = "_")
save(is_heterozygote, file = "./data/heterozygotes_position.Rdata")
#load("./data/heterozygotes_position.Rdata")