source("read_genotypes.R")

library(qtl2)
library(qtl2convert)

pat_snped = inner_join(full_data_Strain, IDs, by = "ID") %>%
  dplyr::select(ID, Strain, Litter_ID_new:Mat_ID)

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID, contains("growth"), Final_weight)

pat_gen = gen %>%
  dplyr::select(ID, chr, pos, gpos, pat_snped$ID)
ggplot(pat_gen, aes(x = gpos, y = chr)) + geom_point(size = 0.3, alpha = 0.5)

current_line = line_order[2]
getFounders = function(current_line){
  current_line_gen = gen %>% dplyr::select(ID, chr, pos, gpos, as.character(filter(pat_snped, Strain == current_line)$ID))
  qtl2convert::find_consensus_geno(as.matrix(current_line_gen[,-c(1:4)]), na.strings = "./.")
}
founders = gsub("/", "", laply(line_order, getFounders))
dimnames(founders)[[1]] = line_order
dimnames(founders)[[2]] = gen$ID
code = rbind(c("1", "0"),
             c("1", "0"),
             c("1", "0"),
             c("1", "0"),
             c("1", "0"),
             c("1", "0"))

founders[,1]
founders_df[1,]
founders_df = data.frame(id = gen$ID, t(encode_geno(founders, code)))
#write2csv(founders_df, "./data/qtl2/founders.csv", overwrite = TRUE)

f6_gen = gen %>% dplyr::select(ID, chr, pos, gpos, as.character(f6_snped$ID))

f6_df = f6_gen[,-c(2:4)] %>%
  replace(., . == "./.", "-") %>%
  replace(., . == "0/0", "B") %>%
  replace(., . == "1/1", "A") %>%
  replace(., . == "0/1", "H")
write2csv(f6_df, "./data/qtl2/f6.csv")


write2csv(select(f6_snped, ID, Final_weight, contains("growth")), "./data/qtl2/pheno.csv")

write2csv(select(f6_snped, ID, Sex) %>% mutate(ngen = 6), "./data/qtl2/covar.csv", overwrite = TRUE)

write2csv(select(gen, ID, chr, gpos), "./data/qtl2/gmap.csv")
write2csv(select(gen, ID, chr, pos), "./data/qtl2/pmap.csv")

f6_qtl2 <- read_cross2("./data/qtl2/control_file.yaml")
