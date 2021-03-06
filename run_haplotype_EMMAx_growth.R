source("read_phenotypes.R")
source("read_genotypes.R")
source("data/LASSO_cluster/emma.r")

if(!require(foreach)){install.packages("foreach"); library(foreach)}
if(!require(doMC)){install.packages("doMC"); library(doMC)}
registerDoMC(2)

load("./data/f6_model_matrices.Rdata")

f6_snped_growth = f6_snped[match(dimnames(f6_model_matrices)[[1]], f6_snped$ID),]
all(f6_snped_growth$ID == dimnames(f6_model_matrices)[[1]])


Xo <- cbind(rep(1, nrow(f6_snped_growth)),.n(f6_snped_growth$Sex)-1)
#K_norm<-(n-1)/sum((diag(n)-matrix(1,n,n)/n)*Af6)*Af6
ids = as.character(f6_snped_growth$ID)
Af6 = (A[ids,ids])
n = nrow(Af6)
Af6_snp = A_snp[f6_snped_growth$ID, f6_snped_growth$ID]
Af6_loco = A_loco[,f6_snped_growth$ID, f6_snped_growth$ID]

source("hapGWAS.R")

growth_traits = colnames(select(f6_snped_growth, contains("growth")))
n_traits = length(growth_traits)

gwas_emmax_ped = vector("list", length(growth_traits))
gwas_emmax_snp = vector("list", length(growth_traits))
gwas_emmax_loco = vector("list", length(growth_traits))
k = 1
for(i in k:n_traits){
  print(growth_traits[[i]])
  gwas_emmax_ped[[i]] = hapGWAS(Y = as.data.frame(f6_snped_growth)[growth_traits[[i]]][,1], K_norm = Af6)
  gwas_emmax_snp[[i]] = hapGWAS(Y = as.data.frame(f6_snped_growth)[growth_traits[[i]]][,1], K_norm = Af6_snp)
  gwas_emmax_loco[[i]] = hapGWASLOCO(Y = as.data.frame(f6_snped_growth)[growth_traits[[i]]][,1], K = Af6_loco)
  k = k + 1
}


save(gwas_emmax_snp, gwas_emmax_ped, gwas_emmax_loco, file = "./data/haplotype_gwas_growth.Rdata")
gwas_emmax_plots = vector("list", length(growth_traits))
for (i in 1:n_traits){
  gwas_bwt_p_ped = ggman(gwas_emmax_ped[[i]], snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt",
                         relative.positions = TRUE, title = "Hap EMMAx PED", sigLine = -log10(2.6e-5),
                         pointSize = 1) + scale_y_continuous(limits = c(0, 13))
  gwas_bwt_p_snp = ggman(gwas_emmax_snp[[i]], snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt",
                         relative.positions = TRUE, title = "Hap EMMAx SNP", sigLine = -log10(2.6e-5),
                         pointSize = 1) + scale_y_continuous(limits = c(0, 13))
  gwas_bwt_p_loco = ggman(gwas_emmax_loco[[i]], snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt",
                          relative.positions = TRUE, title = "Hap EMMAx SNP LOCO", sigLine = -log10(2.6e-5),
                          pointSize = 1) + scale_y_continuous(limits = c(0, 13))

  gwas_emmax_plots[[i]] = plot_grid(gwas_bwt_p_ped, gwas_bwt_p_snp, gwas_bwt_p_loco, ncol = 1)
}
gwas_emmax_plots[[8]]
