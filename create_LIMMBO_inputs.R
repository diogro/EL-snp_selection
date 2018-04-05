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
  replace(., . == "BB", "2") %>%
  replace(., . == "AB", "1")

for(i in 1:dim(df)[1]){
  missing = which(is.na(df[i,]))
  if(length(missing) > 0)
    for(j in 1:length(missing)){
      ind = colnames(df)[missing[j]]
      litter = f6_snped$Litter_ID_new[f6_snped$pID == ind]
      df[i, missing[j]] = round(mean(as.numeric(df[missing[j],filter(f6_snped, Litter_ID_new == litter)$pID]), na.rm = TRUE), 0)
    }
}

phenotypes = na.omit(select(f6_snped, ID, pID, Sex, Final_weight, contains("growth")))

mean_genotype = df %>%
  mutate(snp = paste(chr, pos, paste0("rs", pos), sep = "-")) %>%
  select(snp, phenotypes$pID)


write_csv(mean_genotype, "./data/limmbo_mean_genotype.csv", col_names = TRUE)
write_csv(select(phenotypes, pID, Final_weight, growth_D14D21), "./data/limmbo_phenotypes.csv", col_names = TRUE)

A = 2*kinship(pedAll)
ids = as.character(phenotypes$ID)
Af6 = tbl_df(A[ids,ids])
colnames(Af6) = rownames(Af6) = phenotypes$pID
bend_Af6 = nearPD(as.matrix(Af6))
write_csv(tbl_df(as.matrix(bend_Af6$mat)), "./data/limmbo_relatedness.csv", col_names = TRUE)

write_csv(mutate(select(phenotypes, pID, Sex), Sex = .n(Sex)), "./data/limmbo_covariates.csv", col_names = TRUE)

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

