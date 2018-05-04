source("read_phenotypes.R")

if(!require(vcfR)){install.packages("vcfR"); library(vcfR)}
vcf <- read.vcfR("./data/plink_files/imputed.vcf")