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

fam_file = read_delim("./data/plink_files/per_chrom/atchely_imputed_chr1.fam", delim = " ",
                      col_names = c("litter", "ID",  "sire", "dam", "sex", "pheno"))

fam_pheno = inner_join(select(fam_file, -pheno),
                        select(phenotypes, ID, Final_weight) %>% distinct(ID, .keep_all = TRUE))

write_tsv(select(fam_pheno, litter, ID), "./data/gemma/keep_indviduals.txt", col_names = FALSE)
for(i in 1:20){
  system(paste0("plink --bfile ./data/plink_files/per_chrom/atchely_imputed_chr", i, " --keep ./data/gemma/keep_indviduals.txt --make-bed --out ./data/gemma/growth_chr", i))
  system(paste0("plink --bfile ./data/plink_files/per_chrom/atchely_imputed_not_chr", i, " --keep ./data/gemma/keep_indviduals.txt --make-bed --out ./data/gemma/growth_not_chr", i))
}

fam_file = read_delim("./data/gemma/growth_chr1.fam", delim = " ",
                      col_names = c("litter", "ID",  "sire", "dam", "sex", "pheno"))
fam_pheno = inner_join(select(fam_file, -pheno),
                       select(phenotypes, ID,
                              Final_weight,
                              Litter_size_birth,
                              Birth_litter_size_weaning,
                              Foster_litter_size_weaning) %>% distinct(ID, .keep_all = TRUE))
write_tsv(data.frame(1, fam_pheno$sex,
                     fam_pheno$Litter_size_birth,
                     fam_pheno$Birth_litter_size_weaning,
                     fam_pheno$Foster_litter_size_weaning), "./data/gemma/gemma_covariates.tsv", col_names = FALSE)
for(i in 1:20){
  write_delim(select(fam_pheno, litter:Final_weight), paste0("./data/gemma/growth_chr", i, ".fam"), col_names = FALSE, delim = " ")
  write_delim(select(fam_pheno, litter:Final_weight), paste0("./data/gemma/growth_not_chr", i, ".fam"), col_names = FALSE, delim = " ")
}


for(i in 1:20){
  system(paste0("gemma -bfile data/gemma/growth_not_chr", i, " -gk 1 -o gemma_relatedness_chr", i))
  rel_mat = as.matrix(read_delim(paste0("output/gemma_relatedness_chr", i, ".cXX.txt"), delim = "\t", col_names = FALSE))
  diag(rel_mat) = diag(rel_mat) + 1e-4
  write_delim(x = tbl_df(rel_mat), paste0("output/gemma_relatedness_chr", i, ".cXX.txt"), col_names = FALSE)
}
system("gemma -bfile data/gemma/growth -gk 1 -o gemma_relatedness_chr")



A = 2*kinship(pedigree)
ids = as.character(fam_pheno$ID)
Af6 = (A[ids,ids])
Af6 = Af6 + diag(nrow(Af6)) * 1e-4
colnames(Af6) = rownames(Af6) = phenotypes$ID
bend_Af6 = nearPD(as.matrix(Af6))
colnames(bend_Af6$mat) = rownames(bend_Af6$mat) = phenotypes$ID
write_tsv(tbl_df(as.matrix(bend_Af6$mat)), "./data/gemma/gemma_relatedness.tsv", col_names = FALSE)

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

for(i in 1:20){
system(paste0("gemma \\
        -bfile ./data/gemma/growth_chr", i," \\
        -k output/gemma_relatedness_chr", i, ".cXX.txt \\
        -c ./data/gemma/gemma_covariates.tsv \\
        -lmm 2 \\
        -o growth_r-snp_chr", i))
}

for(i in 1:20){
  system(paste0("gemma \\
        -bfile ./data/gemma/growth_chr", i," \\
        -k data/gemma/gemma_relatedness.tsv \\
        -c ./data/gemma/gemma_covariates.tsv \\
        -lmm 2 \\
        -o growth_r-ped_chr", i))
}

for(i in 1:20){
  system(paste0("gemma \\
        -bfile ./data/gemma/growth_chr", i," \\
        -c ./data/gemma/gemma_covariates.tsv \\
        -lm 2 \\
        -o growth_lm_chr", i))
}

gwas_rsnp = ldply(1:20, function(i) read_tsv(paste0("./output/growth_r-snp_chr",i,".assoc.txt")))
gwas_rped = ldply(1:20, function(i) read_tsv(paste0("./output/growth_r-ped_chr",i,".assoc.txt")))
gwas_lm   = ldply(1:20, function(i) read_tsv(paste0("./output/growth_lm_chr",i,".assoc.txt")))

gwas_qtl_rel = data.frame(rs = names(lrt$p), p_lrt = lrt$p, ps = map$phyPos, chr = map$chr)

qvalue_correction = qvalue(gwas_rsnp$p_lrt, fdr.level = 0.01)
gwas_rsnp = gwas_rsnp %>%
          mutate(qvalues = qvalue_correction$qvalues,
                 significant = qvalue_correction$significant)

qvalue_correction = qvalue(gwas_rped$p_lrt, fdr.level = 0.01)
gwas_rped = gwas_rped %>%
  mutate(qvalues = qvalue_correction$qvalues,
         significant = qvalue_correction$significant)

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

table(gwas$p_lrt < 1e-6)
table(gwas$significant)
hist(gwas_rped$p_lrt)
obs_rsnp      = -log10(sort(gwas_rsnp$p_lrt))
obs_rped      = -log10(sort(gwas_rped$p_lrt))
obs_lm        = -log10(sort(gwas_lm$p_lrt))
obs_qtl_rel   = -log10(sort(lrt$p))
obs_happy     = -log10(sort(gwh$p))
expected = -log10(stats::ppoints(nrow(gwas_rsnp)))
exp_happy =-log10(stats::ppoints(nrow(gwh)))
plot(obs_rped~expected, pch = 19, ylim = c(0, 15))
points(obs_rsnp~expected, pch = 19, col = "blue")
points(obs_qtl_rel~expected, pch = 19, col = "red")
#points(obs_happy~exp_happy, pch = 19, col = "darkgreen")
#points(obs_lm~expected, pch = 19, col = "orange")
abline(0, 1)
abline(h = -log10(5.17E-7))
cov(-log10(gwas_rsnp$p_lrt), -log10(gwas_rped$p_lrt))

hist(gwas_rsnp$p_lrt)
hist(gwas_rped$p_lrt)

table(gwas_rped$p_lrt < 5.17E-7)
gwas_rped[which(gwas_rped$p_lrt < 5.17E-7),]
(gwas_growth_p_ped = ggman(gwas_rped, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt", relative.positions = TRUE, title = "GEMMA ped BWT", sigLine = -log10(2.6e-5), pointSize = 1))
(gwas_growth_p_snp = ggman(gwas_rsnp, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt", relative.positions = TRUE, title = "GEMMA snp BWT", sigLine = -log10(2.6e-5), pointSize = 1))
(gwas_growth_p_qtlRel = ggman(gwas_qtl_rel, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt", relative.positions = TRUE, title = "QTL Rel BWT", sigLine = -log10(2.6e-5), pointSize = 1))
(gwas_growth_p_happy = ggman(gwh, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p", relative.positions = TRUE, title = "Happy BWT", sigLine = -log10(2.6e-5), pointSize = 1))

(gwas_growth_qvalue = ggman(gwas, snp = "rs", bp = "ps", chrom = "chr", pvalue = "qvalues", relative.positions = TRUE, title = "Growth 3 intervals"))

save_plot("~/Dropbox/labbio/data/Atchley project/Genotypes/final_weight_gwas.png", gwas_growth, base_height = 6, base_aspect_ratio = 2)

