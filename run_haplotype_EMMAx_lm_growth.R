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

Y = as.data.frame(f6_snped_growth)[growth_traits[[1]]][,1]
K_norm = Af6_snp

hapGWAS <- function(Y, K_norm) {
  if(any(is.na(Y))) Y[is.na(Y)] = mean(Y, na.rm = TRUE)
  null<-emma.REMLE(Y, as.matrix(Xo), K_norm)

  #EMMAX

  M <- solve(chol(null$vg * K_norm + null$ve * diag(n)))
  Y_t <- crossprod(M, Y)
  Xo_t <- crossprod(M, Xo)

  H0 <- lm(Y_t ~ Xo_t - 1)
  pval = vector("numeric", dim(f6_model_matrices)[[3]])
  pb = txtProgressBar(min = 1, max = length(pval), style = 3)
  k = 1
  for(i in k:length(pval)){
    setTxtProgressBar(pb, i)
    X = f6_model_matrices[,,i]
    X_t <- crossprod(M,X)
    model = lm(Y_t ~ cbind(Xo_t, X_t) - 1)
    comp = anova(H0, model)
    pval[i] = comp$`Pr(>F)`[2]
  }
  close(pb)
  data.frame(rs = gen$ID, ps = gen$pos, chr = gen$chr, p_lrt = pval)
}
hapGWASLOCO <- function(Y, K) {
  if(any(is.na(Y))) Y[is.na(Y)] = mean(Y, na.rm = TRUE)

  #EMMAX
  null = list()
  Y_t = array(0, dim = c(length(Y), 20))
  Xo_t = array(0, dim = c(20, dim(Xo)))
  M = array(0, dim = dim(K))
  H0 = vector("list", 20)

  for(chr in 1:20){
    null[[chr]] <- emma.REMLE(Y, as.matrix(Xo), K[chr,,])
    M[chr,,] <- solve(chol(null[[chr]]$vg*K[chr,,]+null[[chr]]$ve*diag(n)))
    Y_t[,chr] <- crossprod(M[chr,,],Y)
    Xo_t[chr,,] <- crossprod(M[chr,,], Xo)
    H0[[chr]] <- lm(Y_t[,chr] ~ Xo_t[chr,,] - 1)
  }

  pval = vector("numeric", dim(f6_model_matrices)[[3]])
  pb = txtProgressBar(min = 1, max = length(pval), style = 3)
  k = 1
  for(i in k:length(pval)){
    setTxtProgressBar(pb, i)
    chr = gen$chr[i]
    X = f6_model_matrices[,,i]
    X_t <- crossprod(M[chr,,], X)
    model = lm(Y_t[,chr] ~ cbind(Xo_t[chr,,], X_t) - 1)
    comp = anova(H0[[chr]], model)
    pval[i] = comp$`Pr(>F)`[2]
  }
  close(pb)
  data.frame(rs = gen$ID, ps = gen$pos, chr = gen$chr, p_lrt = pval)
}

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
load("./data/haplotype_gwas_growth.Rdata")

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
gwas_emmax_plots[[1]]
positions$chr2 = NULL
positions = dplyr::rename(positions, ps = pos)

genome_wide_p = 2.6e-5

gwas = gwas_emmax_loco[[1]]
gwas = inner_join(gwas, positions, by = c("chr","ps")) %>%
  filter(p_lrt < genome_wide_p)

getChoosenSNPs <- function(gwas) {
  gwas = inner_join(gwas, positions, by = c("chr","ps")) %>%
    filter(p_lrt < genome_wide_p)
  choosen_snp_list = list()
  k = 1
  for(c in unique(gwas$chr)){
    while(TRUE){
      min_p = which.min(gwas[gwas$chr == c, "p_lrt"])
      choosen_snp = as.character(gwas[gwas$chr == c,"rs"][min_p])
      if(length(choosen_snp) > 0){
        choosen_snp_list[[k]] = choosen_snp
        k = k + 1
      } else {
        break
      }
      choosen_pos = gwas[gwas$rs == choosen_snp, "gpos"]
      removable_snps = gwas[gwas$chr == c, "rs"][abs(gwas[gwas$chr == c, "gpos"] - choosen_pos) < 20]
      gwas = filter(gwas, !rs %in% removable_snps)
    }
  }
  unlist(choosen_snp_list)
}
gwas_emmax_ped[[4]] %>% filter(rs %in% getChoosenSNPs(.))
x %>% filter(rs %in% getChoosenSNPs(.))
plot_grid(ggman(gwas_emmax_loco[[1]], snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt",
      relative.positions = TRUE, title = "Hap EMMAx PED", sigLine = -log10(1e-4),
      pointSize = 1) + scale_y_continuous(limits = c(0, 13)),
ggman(x, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt",
      relative.positions = TRUE, title = "Hap EMMAx PED", sigLine = -log10(1e-4),
      pointSize = 1) + scale_y_continuous(limits = c(0, 13)), ncol= 1)

names(gwas_emmax_ped[[1]])
gwas_emmax_loco[[3]] %>%
  filter(p_lrt < 2.6e-5)
