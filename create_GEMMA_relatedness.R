source("read_genotypes.R")

for(i in 1:20){
    fam_file = read_delim(paste0("./data/plink_files/per_chrom/atchley_imputed_not_chr", i, ".fam"), delim = " ",
                          col_names = c("litter", "ID",  "sire", "dam", "sex", "pheno"))
    fam_file$pheno = 1
    write_delim(fam_file, paste0("./data/plink_files/per_chrom/atchley_imputed_not_chr", i, ".fam"), col_names = FALSE, delim = " ")
}
fam_file = read_delim("./data/plink_files/atchley_imputed.fam", delim = " ",
                      col_names = c("litter", "ID",  "sire", "dam", "sex", "pheno"))
fam_file$pheno = 1
write_delim(fam_file, "./data/plink_files/atchley_imputed.fam", col_names = FALSE, delim = " ")

registerDoMC(20)
foreach(i=1:20) %dopar% {
  system(paste0("gemma -bfile ./data/plink_files/per_chrom/atchley_imputed_not_chr", i, " -gk 1 -o gemma_relatedness_full_chr", i))
}
system("gemma -bfile ./data/plink_files/atchley_imputed -gk 1 -o gemma_relatedness_full")

n = nrow(fam_file)
A_loco = array(0, dim = c(20, n, n))
for(i in 1:20){
    A_loco[i,,] = as.matrix(read_tsv(paste0("./output/gemma_relatedness_full_chr", i, ".cXX.txt"), col_names = FALSE)) + diag(n) * 1e-4
}
A_snp = as.matrix(read_tsv("./output/gemma_relatedness_full.cXX.txt", col_names = FALSE)) + diag(n) * 1e-4
fam_file = read_delim(paste0("./data/plink_files/per_chrom/atchley_imputed_not_chr", i, ".fam"), delim = " ",
                      col_names = c("litter", "ID",  "sire", "dam", "sex", "pheno"))
dimnames(A_loco) = list(1:20, fam_file$ID, fam_file$ID)
dimnames(A_snp) = list(fam_file$ID, fam_file$ID)
save(A_snp, A_loco, file = "./data/gemma_relatedness.Rdata")
