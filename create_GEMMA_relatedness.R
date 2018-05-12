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
