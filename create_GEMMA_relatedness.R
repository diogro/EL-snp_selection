source("read_genotypes.R")

registerDoMC(20)
foreach(i=1:20) %dopar% {
  system(paste0("gemma -bfile ./data/plink_files/per_chrom/atchley_imputed_not_chr", i, " -gk 1 -o gemma_relatedness_full_chr", i))

}
system("gemma -bfile ./data/plink_files/atchley_imputed -gk 1 -o gemma_relatedness_full")