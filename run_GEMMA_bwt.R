source("read_genotypes.R")

f5f6_snped = inner_join(full_data_F5F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Foster_litter_size_weaning, contains("growth"), Final_weight)

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Foster_litter_size_weaning, contains("growth"), Final_weight)

phenotypes = select(f6_snped, ID, Sex, Final_weight,
                    Litter_size_birth,
                    Birth_litter_size_weaning,
                    Foster_litter_size_weaning,
                    growth_D3D7:growth_D21D28) %>% distinct(ID, .keep_all = TRUE)

fam_file = read_delim("./data/plink_files/per_chrom/atchley_imputed_chr1.fam", delim = " ",
                      col_names = c("litter", "ID",  "sire", "dam", "sex", "pheno"))
fam_file$ID = as.character(fam_file$ID)

fam_pheno = inner_join(select(fam_file, -pheno),
                        select(phenotypes, ID, Final_weight) %>% distinct(ID, .keep_all = TRUE))
fam_file$ID = as.numeric(fam_file$ID)

write_tsv(select(fam_pheno, litter, ID), "./data/gemma/keep_indviduals.txt", col_names = FALSE)
system("plink --bfile ./data/plink_files/atchley_imputed --keep ./data/gemma/keep_indviduals.txt --make-bed --out ./data/gemma/bwt")
for(i in 1:20){
  system(paste0("plink --bfile ./data/plink_files/per_chrom/atchley_imputed_chr", i, " --keep ./data/gemma/keep_indviduals.txt --make-bed --out ./data/gemma/bwt_chr", i))
  system(paste0("plink --bfile ./data/plink_files/per_chrom/atchley_imputed_not_chr", i, " --keep ./data/gemma/keep_indviduals.txt --make-bed --out ./data/gemma/bwt_not_chr", i))
}

fam_file = read_delim("./data/gemma/bwt_chr1.fam", delim = " ",
                      col_names = c("litter", "ID",  "sire", "dam", "sex", "pheno"))
fam_file$ID = as.character(fam_file$ID)
fam_pheno = inner_join(select(fam_file, -pheno),
                       select(phenotypes, ID,
                              Final_weight,
                              Litter_size_birth,
                              Birth_litter_size_weaning,
                              Foster_litter_size_weaning) %>% distinct(ID, .keep_all = TRUE))
fam_pheno$ID = as.numeric(fam_pheno$ID)
write_tsv(data.frame(1, fam_pheno$sex,
                     fam_pheno$Litter_size_birth,
                     fam_pheno$Birth_litter_size_weaning,
                     fam_pheno$Foster_litter_size_weaning), "./data/gemma/gemma_covariates.tsv", col_names = FALSE)
for(i in 1:20){
  write_delim(select(fam_pheno, litter:Final_weight), paste0("./data/gemma/bwt_chr", i, ".fam"), col_names = FALSE, delim = " ")
  write_delim(select(fam_pheno, litter:Final_weight), paste0("./data/gemma/bwt_not_chr", i, ".fam"), col_names = FALSE, delim = " ")
}
write_delim(select(fam_pheno, litter:Final_weight), "./data/gemma/bwt.fam", col_names = FALSE, delim = " ")

registerDoMC(8)
for(i in 1:20)  {
  system(paste0("gemma -bfile data/gemma/bwt_not_chr", i, " -gk 1 -o gemma_relatedness_chr", i))
  #rel_mat = as.matrix(read_delim(paste0("output/gemma_relatedness_chr", i, ".cxx.txt"), delim = "\t", col_names = false))
  #diag(rel_mat) = diag(rel_mat) + 1e-4
  #write_delim(x = tbl_df(rel_mat), paste0("output/gemma_relatedness_chr", i, ".cXX.txt"), col_names = FALSE)
}
system("gemma -bfile data/gemma/bwt -gk 1 -o gemma_relatedness")

A = 2*kinship(pedigree)
ids = as.character(fam_pheno$ID)
Af6 = (A[ids,ids])
Af6 = Af6 + diag(nrow(Af6)) * 1e-4
colnames(Af6) = rownames(Af6) = fam_pheno$ID
write_tsv(tbl_df(as.matrix(Af6)), "./data/gemma/gemma_relatedness.tsv", col_names = FALSE)

# phenotypes$animal = phenotypes$ID
# Ainv = as(base::solve(bend_Af6$mat), "dgCMatrix")
# colnames(Ainv) = rownames(Ainv) = phenotypes$ID
# prior_multi <- list(G = list(G1 = list(V = diag(3), nu = 4)),
#                   R = list(V = diag(3), nu = 0.002))
# null_mcmc = MCMCglmm(cbind(bwt_D3D7, bwt_D7D14, bwt_D21D28) ~ trait:Sex - 1,
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

library(foreach)
library(doMC)
registerDoMC(20)
foreach(i=1:20) %dopar% {
system(paste0("gemma \\
        -bfile ./data/gemma/bwt_chr", i," \\
        -k output/gemma_relatedness_chr", i, ".cXX.txt \\
        -c ./data/gemma/gemma_covariates.tsv \\
        -lmm 2 \\
        -o bwt_r-LOCO_snp_chr", i))
}

foreach(i=1:20) %dopar% {
system(paste0("gemma \\
        -bfile ./data/gemma/bwt_chr", i," \\
        -k output/gemma_relatedness.cXX.txt \\
        -c ./data/gemma/gemma_covariates.tsv \\
        -lmm 2 \\
        -o bwt_r-snp_chr", i))
}

foreach(i=1:20) %dopar% {
  system(paste0("gemma \\
        -bfile ./data/gemma/bwt_chr", i," \\
        -k data/gemma/gemma_relatedness.tsv \\
        -c ./data/gemma/gemma_covariates.tsv \\
        -lmm 2 \\
        -o bwt_r-ped_chr", i))
}

gwas_rsnp = ldply(1:20, function(i) read_tsv(paste0("data/output/bwt_r-snp_chr",i,".assoc.txt")))
gwas_rsnp_loco = ldply(1:20, function(i) read_tsv(paste0("data/output/bwt_r-LOCO_snp_chr",i,".assoc.txt")))
gwas_rped = ldply(1:20, function(i) read_tsv(paste0("data/output/bwt_r-ped_chr",i,".assoc.txt")))


table(gwas_rped$p_lrt < 5.17E-7)
gwas_rped[which(gwas_rped$p_lrt < 5.17E-7),]
gwas_bwt_GEMMA_ped = ggman(gwas_rped, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt", relative.positions = TRUE, title = "GEMMA bwt \n Pedigree relatedness", sigLine = -log10(2.6e-5), pointSize = 0.5) + scale_y_continuous(limits = c(0, 13))
gwas_bwt_GEMMA_loco = ggman(gwas_rsnp_loco, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt", relative.positions = TRUE, title = "GEMMA bwt \n SNP relatedness LOCO", sigLine = -log10(2.6e-5), pointSize = 0.5) + scale_y_continuous(limits = c(0, 13))
gwas_bwt_GEMMA_snp = ggman(gwas_rsnp, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt", relative.positions = TRUE, title = "GEMMA bwt \n SNP relatedness", sigLine = -log10(2.6e-5), pointSize = 0.5) + scale_y_continuous(limits = c(0, 13))
all_gemma = plot_grid(gwas_bwt_GEMMA_ped, gwas_bwt_GEMMA_snp, gwas_bwt_GEMMA_loco, ncol = 3)

save_plot("~/gemma_bwt_ped.png", gwas_bwt_GEMMA_ped, base_height = 7, base_aspect_ratio = 2)
save_plot("~/gemma_bwt_loco.png", gwas_bwt_GEMMA_loco, base_height = 7, base_aspect_ratio = 2)
save_plot("~/gemma_bwt_snp.png", gwas_bwt_GEMMA_snp, base_height = 7, base_aspect_ratio = 2)
save_plot("~/gemma_bwt_all.png", all_gemma, base_height = 6, base_aspect_ratio = 1, ncol = 3, nrow = 1)
