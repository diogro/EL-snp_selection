source("read_genotypes.R")

f5f6_snped = inner_join(full_data_F5F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Foster_litter_size_weaning, contains("growth"), Final_weight)

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Foster_litter_size_weaning, contains("growth"), Final_weight)

phenotypes = select(f6_snped, ID, Sex, Final_weight,
                    Litter_size_birth,
                    Birth_litter_size_weaning,
                    Foster_litter_size_weaning,
                    growth_traits) %>% distinct(ID, .keep_all = TRUE)
phenotypes[,growth_traits] = scale(phenotypes[,growth_traits])

fam_file = read_delim("./data/plink_files/per_chrom/atchley_imputed_chr1.fam", delim = " ",
                      col_names = c("litter", "ID",  "sire", "dam", "sex", "pheno"))

fam_file$ID = as.character(fam_file$ID)
fam_pheno = inner_join(select(fam_file, -pheno),
                        select(phenotypes, ID, contains("growth")) %>% distinct(ID, .keep_all = TRUE))
fam_pheno$ID = as.numeric(fam_pheno$ID)

write_tsv(select(fam_pheno, litter, ID), "./data/gemma/keep_indviduals.txt", col_names = FALSE)
for(i in 1:20){
  system(paste0("plink --bfile ./data/plink_files/per_chrom/atchley_imputed_chr", i, " --keep ./data/gemma/keep_indviduals.txt --make-bed --out ./data/gemma/growth_chr", i))
  system(paste0("plink --bfile ./data/plink_files/per_chrom/atchley_imputed_not_chr", i, " --keep ./data/gemma/keep_indviduals.txt --make-bed --out ./data/gemma/growth_not_chr", i))
}

fam_file = read_delim("./data/gemma/growth_chr1.fam", delim = " ",
                      col_names = c("litter", "ID",  "sire", "dam", "sex", "pheno"))
fam_file$ID = as.character(fam_file$ID)
fam_pheno = inner_join(select(fam_file, -pheno),
                       select(phenotypes, ID,
                              contains("growth"),
                              Litter_size_birth,
                              Birth_litter_size_weaning,
                              Foster_litter_size_weaning) %>% distinct(ID, .keep_all = TRUE))
fam_pheno$ID = as.numeric(fam_pheno$ID)
write_tsv(data.frame(1, fam_pheno$sex,
                     fam_pheno$Litter_size_birth,
                     fam_pheno$Birth_litter_size_weaning,
                     fam_pheno$Foster_litter_size_weaning), "./data/gemma/gemma_covariates.tsv", col_names = FALSE)
for(i in 1:20){
  write_delim(select(fam_pheno, litter:growth_D42D56), paste0("./data/gemma/growth_chr", i, ".fam"), col_names = FALSE, delim = " ")
  write_delim(select(fam_pheno, litter:growth_D42D56), paste0("./data/gemma/growth_not_chr", i, ".fam"), col_names = FALSE, delim = " ")
}



library(foreach)
library(doMC)
registerDoMC(20)
foreach(i=1:20) %dopar% {
 system(paste0("gemma -bfile data/gemma/growth_not_chr", i, " -gk 1 -o gemma_relatedness_chr", i))
## rel_mat = as.matrix(read_delim(paste0("output/gemma_relatedness_chr", i, ".cXX.txt"), delim = "\t", col_names = FALSE))
## diag(rel_mat) = diag(rel_mat) + 1e-3
#  write_delim(x = tbl_df(rel_mat), paste0("output/gemma_relatedness_chr", i, ".cXX.txt"), col_names = FALSE)
}
#system("gemma -bfile data/gemma/growth -gk 1 -ogemma_relatedness_chr")

A = 2*kinship(pedigree)
ids = as.character(fam_pheno$ID)
Af6 = (A[ids,ids])
Af6 = Af6 + diag(nrow(Af6)) * 1e-4
colnames(Af6) = rownames(Af6) = phenotypes$ID
write_tsv(tbl_df(Af6), "./data/gemma/pedigree_relatedness.tsv", col_names = FALSE)

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

registerDoMC(20)
foreach(i=1:20) %dopar% {
system(paste0("gemma \\
        -bfile ./data/gemma/growth_chr", i," \\
        -k data/gemma/pedigree_relatedness.tsv \\
        -c ./data/gemma/gemma_covariates.tsv \\
        -lmm 4 \\
        -n 1 2 3 4\\
        -o growth_r-ped_multivariate_2w_t1-4_chr", i))
  sprintf("SNP PED %d, traits 1-4", i)
}

registerDoMC(20)
foreach(i=1:20) %dopar% {
system(paste0("gemma \\
        -bfile ./data/gemma/growth_chr", i," \\
        -k output/gemma_relatedness_chr", i, ".cXX.txt \\
        -c ./data/gemma/gemma_covariates.tsv \\
        -lmm 4 \\
        -n 1 2 3 4\\
        -o growth_r-loco_multivariate_2w_t1-4_chr", i))
  sprintf("SNP loco %d, traits 1-4", i)
}

