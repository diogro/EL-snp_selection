source("read_genotypes.R")

f5f6_snped = inner_join(full_data_F5F6, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID, contains("growth"), Final_weight)

f5_snped = inner_join(full_data_F5, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)

f6_genotypes = gen %>% select(ID, chr, pos, gpos, f6_snped$pID)

df_gemma = f6_genotypes %>%
  replace(., . == "NoCall", NA) %>%
  replace(., . == "AA", "0") %>%
  replace(., . == "BB", "1") %>%
  replace(., . == "AB", "0.5")

for(i in 1:dim(df_gemma)[1]){
  missing = which(is.na(df_gemma[i,]))
  if(length(missing) > 0)
    for(j in 1:length(missing)){
      ind = colnames(df_gemma)[missing[j]]
      litter = f6_snped$Litter_ID_new[f6_snped$pID == ind]
      df_gemma[i, missing[j]] = mean(as.numeric(df_gemma[missing[j],filter(f6_snped, Litter_ID_new == litter)$pID]), na.rm = TRUE)
    }
}

phenotypes = select(f6_snped, ID, pID, Sex, growth_D3D7:growth_D14D21)

mean_genotype = df_gemma %>%
  mutate(chr = "B", pos = "A") %>%
  rename(minor = chr,
         major = pos) %>%
  select(ID, minor, major, -gpos, f6_snped$pID)
write_csv(mean_genotype, "./data/gemma_mean_genotype.csv", col_names = FALSE)
write_tsv(select(phenotypes, growth_D3D7:growth_D14D21), "./data/gemma_phenotypes_growth_D2-D7-D21.tsv", col_names = FALSE)
write_csv(select(df_gemma, ID, pos, chr), "./data/gemma_marker_position.csv", col_names = FALSE)

A = 2*kinship(pedAll)
ids = as.character(phenotypes$ID)
Af6 = tbl_df(A[ids,ids])
colnames(Af6) = rownames(Af6) = phenotypes$pID
bend_Af6 = nearPD(as.matrix(Af6))
write_tsv(tbl_df(as.matrix(bend_Af6$mat)), "./data/gemma_relatedness.tsv", col_names = FALSE)

write_tsv(data.frame(1, .n(phenotypes$Sex)), "./data/gemma_covariates.tsv", col_names = FALSE)

#phenotypes$animal = phenotypes$ID
#ped2$animal = ped2$ID
#prior_multi <- list(G = list(G1 = list(V = diag(3), nu = 4)),
                  #R = list(V = diag(3), nu = 0.002))
#null_mcmc = MCMCglmm(cbind(growth_D3D7, growth_D7D14, growth_D14D21) ~ trait:Sex - 1,
                     #random = ~us(trait):animal,
                     #rcov = ~us(trait):units,
                     #prior = prior_multi,
                     #family = rep("gaussian", 3),
                     #pedigree = ped2[,c("animal", "dam", "sire")],
                     #data = as.data.frame(phenotypes))
#summary(null_mcmc)
#plot(null_mcmc)
#nm = colnames(null_mcmc$VCV)
#G = matrix(colMeans(null_mcmc$VCV)[grep("animal", nm)], 3, 3)
#R = matrix(colMeans(null_mcmc$VCV)[grep("units", nm)], 3, 3)
#cov(phenotypes[,4:6])
#G + R

system("./gemma/gemma.linux \\
        -g ./data/gemma_mean_genotype.csv \\
        -p ./data/gemma_phenotypes_growth_D2-D7-D21.tsv \\
        -a ./data/gemma_marker_position.csv \\
        -k ./data/gemma_relatedness.tsv \\
        -c ./data/gemma_covariates.tsv \\
        -n 1 2 3 \\
        -lmm 4 -o growth")

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

