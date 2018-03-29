source("read_genotypes.R")

f5f6_snped = inner_join(full_data_F5F6, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID, contains("growth"), Final_weight)

f5_snped = inner_join(full_data_F5, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)

f6_genotypes = gen %>% select(ID, chr, pos, gpos, f6_snped$pID)

df = f6_genotypes %>%
  replace(., . == "NoCall", NA) %>%
  replace(., . == "AA", "0") %>%
  replace(., . == "BB", "1") %>%
  replace(., . == "AB", "0.5")

for(i in 1:dim(df)[1]){
  missing = which(is.na(df[i,]))
  if(length(missing) > 0)
    for(j in 1:length(missing)){
      ind = colnames(df)[missing[j]]
      litter = f6_snped$Litter_ID_new[f6_snped$pID == ind]
      df[i, missing] = mean(as.numeric(df[missing[j],filter(f6_snped, Litter_ID_new == litter)$pID]), na.rm = TRUE)
    }
}

mean_genotype = df %>%
  mutate(chr = "B", pos = "A") %>%
  rename(minor = chr,
         major = pos) %>%
  select(ID, minor, major, -gpos, f6_snped$pID)
write_csv(mean_genotype, "./data/gemma_mean_genotype.csv", col_names = FALSE)
write_tsv(select(f6_snped, growth_D0D3:growth_D35D42), "./data/gemma_phenotypes.tsv", col_names = FALSE)
write_csv(select(df, ID, pos, chr), "./data/gemma_marker_position.csv", col_names = FALSE)

A = 2*kinship(pedAll)
ids = as.character(f6_snped$ID)
Af6 = tbl_df(A[ids,ids])
write_tsv(Af6, "./data/gemma_relatedness.tsv", col_names = FALSE)

eigA = eigen(Af6)
plot(eigA$values)

write_tsv(data.frame(1, .n(f6_snped$Sex)), "./data/gemma_covariates.tsv", col_names = FALSE)

library(devtools)
install_github("drveera/ggman")
library(ggman)
if(!require(qvalue)){install.packages("qvalue"); library(qvalue)}

gwas = read_tsv("./gemma/output/output.csv.assoc.txt")

qvalue_correction = qvalue(gwas$p_wald, fdr.level = 0.05)
gwas = gwas %>%
          mutate(qvalues = qvalue_correction$qvalues,
                 significant = qvalue_correction$significant)

gwas_final_weight = ggman(gwas, snp = "rs", bp = "ps", chrom = "chr", pvalue = "qvalues", relative.positions = TRUE, title = "Final Weight", sigLine = 1.3)


save_plot("~/Dropbox/labbio/data/Atchley project/Genotypes/final_weight_gwas.png", gwas_growth, base_height = 6, base_aspect_ratio = 2)

