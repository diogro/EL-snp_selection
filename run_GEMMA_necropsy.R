source("read_genotypes.R")

f5f6_snped = inner_join(full_data_F5F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Foster_litter_size_weaning, contains("growth"), Final_weight)

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Foster_litter_size_weaning, contains("growth"), Final_weight)

phenotypes = select(f6_snped, ID, Sex, Final_weight,
                    Litter_size_birth,
                    Birth_litter_size_weaning,
                    Foster_litter_size_weaning,
                    contains("growth")) %>% distinct(ID, .keep_all = TRUE)

fam_file = read_delim("./data/plink_files/per_chrom/atchely_imputed_chr1.fam", delim = " ",
                      col_names = c("litter", "ID",  "sire", "dam", "sex", "pheno"))

fam_file$ID = as.character(fam_file$ID)
fam_pheno = inner_join(select(fam_file, -pheno),
                        select(phenotypes, ID, contains("growth")) %>% distinct(ID, .keep_all = TRUE))
fam_pheno$ID = as.numeric(fam_pheno$ID)

write_tsv(select(fam_pheno, litter, ID), "./data/gemma/keep_indviduals.txt", col_names = FALSE)
for(i in 1:20){
  system(paste0("plink --bfile ./data/plink_files/per_chrom/atchely_imputed_chr", i, " --keep ./data/gemma/keep_indviduals.txt --make-bed --out ./data/gemma/growth_chr", i))
  system(paste0("plink --bfile ./data/plink_files/per_chrom/atchely_imputed_not_chr", i, " --keep ./data/gemma/keep_indviduals.txt --make-bed --out ./data/gemma/growth_not_chr", i))
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
  write_delim(select(fam_pheno, litter:growth_D49D56), paste0("./data/gemma/growth_chr", i, ".fam"), col_names = FALSE, delim = " ")
  write_delim(select(fam_pheno, litter:growth_D49D56), paste0("./data/gemma/growth_not_chr", i, ".fam"), col_names = FALSE, delim = " ")
}


foreach(i=1:20) %dopar% {
  system(paste0("gemma -bfile data/gemma/growth_not_chr", i, " -gk 1 -o gemma_relatedness_chr", i))
 # rel_mat = as.matrix(read_delim(paste0("output/gemma_relatedness_chr", i, ".cXX.txt"), delim = "\t", col_names = FALSE))
 # diag(rel_mat) = diag(rel_mat) + 1e-3
#  write_delim(x = tbl_df(rel_mat), paste0("output/gemma_relatedness_chr", i, ".cXX.txt"), col_names = FALSE)
}
#system("gemma -bfile data/gemma/growth -gk 1 -o gemma_relatedness_chr")
for(i in 1:20) {
  system(paste0("gemma -bfile data/gemma/growth_not_chr", i, " -gk 1 -o gemma_relatedness_chr", i))
  # rel_mat = as.matrix(read_delim(paste0("output/gemma_relatedness_chr", i, ".cXX.txt"), delim = "\t", col_names = FALSE))
  # diag(rel_mat) = diag(rel_mat) + 1e-3
  #  write_delim(x = tbl_df(rel_mat), paste0("output/gemma_relatedness_chr", i, ".cXX.txt"), col_names = FALSE)
}

A = 2*kinship(pedigree)
ids = as.character(fam_pheno$ID)
Af6 = (A[ids,ids])
Af6 = Af6 + diag(nrow(Af6)) * 1e-3
colnames(Af6) = rownames(Af6) = phenotypes$ID
write_tsv(tbl_df(Af6), "./data/gemma/gemma_relatedness.tsv", col_names = FALSE)

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

library(foreach)
library(doMC)
registerDoMC(10)
foreach(i=1:20) %dopar% {
system(paste0("gemma \\
        -bfile ./data/gemma/growth_chr", i," \\
        -k output/gemma_relatedness_chr", i, ".cXX.txt \\
        -c ./data/gemma/gemma_covariates.tsv \\
        -lmm 2 \\
        -n 1 2 3 \\
        -o growth_r-snp_chr", i))
  sprintf("SNP LOCO %d, traits 1-5", i)
}

gwas_rsnp = ldply(1:20, function(i) read_tsv(paste0("./output/growth_r-snp_chr",i,".assoc.txt")))
gwas_rped = ldply(1:20, function(i) read_tsv(paste0("./output/growth_r-ped_chr",i,".assoc.txt")))

plot.inflation <- function (x, size = 2) {

  # Get the number of p-values.
  n <- length(x)

  # Compute the negative log10(p-values), and sort them from largest
  # to smallest.
  y <- rev(sort(-log10(x)))

  # Create the q-q plot.
  return(ggplot(data.frame(x = -log10((1:n)/n),y = y),aes(x = x,y = y)) +
           geom_abline(intercept = 0,slope = 1,color = "magenta") +
           geom_point(color = "dodgerblue",shape = 20,size = 2) +
           labs(x = "Expected -log10 p-value",y = "Observed -log10 p-value") +
           theme(axis.line = element_blank()))
}

plot.inflation(gwas_rped$p_lrt)
plot.inflation(gwas_rsnp$p_lrt)

table(gwas_rsnp$p_lrt < 5.17E-7)
gwas_rped[which(gwas_rped$p_lrt < 5.17E-7),]
(gwas_growth_p_ped = ggman(gwas_rped, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt", relative.positions = TRUE, title = "GEMMA ped growth", sigLine = -log10(2.6e-5), pointSize = 1))
(gwas_growth_p_snp = ggman(gwas_rsnp, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt", relative.positions = TRUE, title = "GEMMA snp growth", sigLine = -log10(2.6e-5), pointSize = 1))
(gwas_growth_p_qtlRel = ggman(gwas_qtl_rel, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt", relative.positions = TRUE, title = "QTL Rel growth", sigLine = -log10(2.6e-5), pointSize = 1))
(gwas_growth_p_happy = ggman(gwh, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p", relative.positions = TRUE, title = "Happy growth", sigLine = -log10(2.6e-5), pointSize = 1))

(gwas_growth_qvalue = ggman(gwas, snp = "rs", bp = "ps", chrom = "chr", pvalue = "qvalues", relative.positions = TRUE, title = "growth 3 intervals"))

save_plot("~/Dropbox/labbio/data/Atchley project/Genotypes/final_weight_gwas.png", gwas_growth, base_height = 6, base_aspect_ratio = 2)

