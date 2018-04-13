source("read_genotypes.R")

f5f6_snped = inner_join(full_data_F5F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID, contains("growth"), Final_weight)

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID, contains("growth"), Final_weight)

phenotypes = select(f6_snped, ID, Sex, Final_weight, growth_D3D7:growth_D21D28) %>% distinct(ID, .keep_all = TRUE)

fam_file = read_delim("./data/plink_files/atchely_imputed_thinned.fam", delim = " ",
                      col_names = c("litter", "ID",  "sire", "dam", "sex", "pheno"))

fam_pheno = inner_join(select(fam_file, -pheno),
                        select(phenotypes, ID, Final_weight) %>% distinct(ID, .keep_all = TRUE))

write_tsv(select(fam_pheno, litter, ID), "./data/gemma/keep_indviduals.txt", col_names = FALSE)
system("plink --bfile ./data/plink_files/atchely_imputed_thinned --keep ./data/gemma/keep_indviduals.txt --make-bed --out ./data/gemma/growth")

fam_file = read_delim("./data/gemma/growth.fam", delim = " ",
                      col_names = c("litter", "ID",  "sire", "dam", "sex", "pheno"))
fam_pheno = inner_join(select(fam_file, -pheno),
                       select(phenotypes, ID, Final_weight) %>% distinct(ID, .keep_all = TRUE))
write_tsv(data.frame(1, fam_pheno$sex), "./data/gemma/gemma_covariates.tsv", col_names = FALSE)
write_delim(fam_pheno, "./data/gemma/growth.fam", col_names = FALSE, delim = " ")

system("gemma -bfile data/gemma/growth -gk 1 -o gemma_relatedness")

# A = 2*kinship(pedigree)
# ids = as.character(fam_pheno$ID)
# Af6 = (A[ids,ids])
# colnames(Af6) = rownames(Af6) = phenotypes$ID
# bend_Af6 = nearPD(as.matrix(Af6))
# colnames(bend_Af6$mat) = rownames(bend_Af6$mat) = phenotypes$ID
# write_tsv(tbl_df(as.matrix(bend_Af6$mat)), "./data/gemma/gemma_relatedness.tsv", col_names = FALSE)

# phenotypes$animal = phenotypes$ID
# Ainv = as(base::solve(bend_Af6$mat), "dgCMatrix")
# colnames(Ainv) = rownames(Ainv) = phenotypes$ID
# prior_multi <- list(G = list(G1 = list(V = diag(3), nu = 4)),
#                   R = list(V = diag(3), nu = 0.002))
# null_mcmc = MCMCglmm(cbind(growth_D3D7, growth_D7D14, growth_D21D28) ~ trait:Sex - 1,
#                      random = ~us(trait):animal,
#                      rcov = ~us(trait):units,
#                      prior = prior_multi,
#                      family = rep("gaussian", 3),
#                      ginverse = list(animal = Ainv),
#                      data = as.data.frame(phenotypes))
#summary(null_mcmc)
#plot(null_mcmc)
#nm = colnames(null_mcmc$VCV)
#G = matrix(colMeans(null_mcmc$VCV)[grep("animal", nm)], 3, 3)
#R = matrix(colMeans(null_mcmc$VCV)[grep("units", nm)], 3, 3)
#cov(phenotypes[,4:6])
#G + R

system("gemma \\
        -bfile ./data/gemma/growth \\
        -k output/gemma_relatedness.cXX.txt \\
        -c ./data/gemma/gemma_covariates.tsv \\
        -lmm 2 \\
        -o growth")

gwas = read_tsv("./output/growth.assoc.txt")

qvalue_correction = qvalue(gwas$p_lrt, fdr.level = 0.05)
gwas = gwas %>%
          mutate(qvalues = qvalue_correction$qvalues,
                 significant = qvalue_correction$significant)
table(gwas$p_lrt < 1e-3)

(gwas_growth_p_lrt = ggman(gwas, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt", relative.positions = TRUE, title = "Growth 3 intervals", sigLine = 3))
(gwas_growth_qvalue = ggman(gwas, snp = "rs", bp = "ps", chrom = "chr", pvalue = "qvalues", relative.positions = TRUE, title = "Growth 3 intervals"))

save_plot("~/Dropbox/labbio/data/Atchley project/Genotypes/final_weight_gwas.png", gwas_growth, base_height = 6, base_aspect_ratio = 2)

