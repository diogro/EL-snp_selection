source("classify_haplotypes_4_realz.R")
source("data/LASSO_cluster/emma.r")

gen2 = select(gen,  ID, chr, pos, gpos)
rm(gen)
gen = gen2
gc()
Xo <- cbind(.n(f6_snped$Sex)-1,
            scale(cbind(f6_snped$Litter_size_birth,
                  f6_snped$Birth_litter_size_weaning, f6_snped$Foster_litter_size_weaning)))
colnames(Xo) = c("sex", "lsb", "blsw", "flsw")
n = nrow(Af6)

X = f6_model_matrices[,,10]
caret::findLinearCombos(cbind(1, Xo,X[,6:1]))

lm_no_int = lm(scale(Final_weight)~ 0 + cbind(Xo,X), data = f6_snped)
lm_int = lmer(scale(Final_weight)~ cbind(Xo,X) + (1|Mat_ID), data = f6_snped)
c_no_int = summary(lm_no_int)
c_no_int$coefficients[5:10,1]
c_int = summary(lm_int)
c_int$coefficients[6:10,1]

Af6_snp = A_snp[f6_snped$ID, f6_snped$ID]
Af6_loco = A_loco[,f6_snped$ID, f6_snped$ID]

hapGWAS <- function(Y, K) {
  Y[is.na(Y)] = mean(Y, na.rm = TRUE)
  Y = scale(Y)

  #EMMAX

  null<-emma.REMLE(Y,as.matrix(Xo),K)
  M<-solve(chol(null$vg*K+null$ve*diag(n)))
  Y_t<-crossprod(M,Y)
  Xo_t <- crossprod(M,Xo)

  RSSf = vector("numeric", dim(f6_model_matrices)[[3]])
  k = 1
  for(i in k:length(RSSf)){
    X = f6_model_matrices[,,i]
    X_t <- crossprod(M,X)
    RSSf[i] = sum(lsfit(cbind(Xo_t,X_t),Y_t,intercept = FALSE)$residuals^2)
  }
  m = length(RSSf)
  RSS_H0 <- rep(sum(lsfit(Xo_t,Y_t,intercept = FALSE)$residuals^2),m)
  df1<-1
  df2<-n-df1-1
  R2<-1-1/(RSS_H0/RSSf)
  F<-(RSS_H0/RSSf-1)*df2/df1
  pval<-pf(F,df1,df2,lower.tail=FALSE)
  data.frame(rs = gen$ID, ps = gen$pos, chr = gen$chr, p_lrt = pval)
}
hapGWASLOCO <- function(Y, K) {
  Y[is.na(Y)] = mean(Y, na.rm = TRUE)
  Y = scale(Y)

  #EMMAX
  null = list()
  Y_t = array(0, dim = c(length(Y), 20))
  Xo_t = array(0, dim = c(20, dim(Xo)))
  M = array(0, dim = dim(K))
  for(i in 1:20){
    null[[i]] <- emma.REMLE(Y, as.matrix(Xo), K[i,,])
    M[i,,] <- solve(chol(null[[i]]$vg*K[i,,]+null[[i]]$ve*diag(n)))
    Y_t[,i] <- crossprod(M[i,,],Y)
    Xo_t[i,,] <- crossprod(M[i,,], Xo)
  }

  RSSf = vector("numeric", dim(f6_model_matrices)[[3]])
  k = 1
  for(i in k:length(RSSf)){
    chr = gen$chr[i]
    X = f6_model_matrices[,,i]
    X_t <- crossprod(M[chr,,], X)
    RSSf[i] = sum(lsfit(cbind(Xo_t[chr,,],X_t),Y_t[,chr],intercept = FALSE)$residuals^2)
  }
  m = length(RSSf)
  H0 = vector("numeric", 20)
  loci_per_chr = vector("numeric", 20)
  for(i in 1:20){
    H0[i] <- sum(lsfit(Xo_t[i,,],Y_t[,i],intercept = FALSE)$residuals^2)
    loci_per_chr[i] = sum(gen$chr == i)
  }
  RSS_H0 = rep(H0, loci_per_chr)
  df1<-1
  df2<-n-df1-1
  R2<-1-1/(RSS_H0/RSSf)
  F<-(RSS_H0/RSSf-1)*df2/df1
  pval<-pf(F,df1,df2,lower.tail=FALSE)
  data.frame(rs = gen$ID, ps = gen$pos, chr = gen$chr, p_lrt = pval)
}

gwas_emmax_ped = hapGWAS(f6_snped$Final_weight, Af6)
gwas_emmax_snp = hapGWAS(f6_snped$Final_weight, Af6_snp)
gwas_emmax_loco = hapGWASLOCO(f6_snped$Final_weight, Af6_loco)
save(gwas_emmax_snp, gwas_emmax_ped, gwas_emmax_loco, file = "./data/haplotype_gwas_bwt.Rdata")
load("./data/haplotype_gwas_bwt.Rdata")
gwas_bwt_HAP_ped = ggman(gwas_emmax_ped, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt",
                       relative.positions = TRUE, title = "Hap EMMAx PED", sigLine = -log10(2.6e-5),
                       pointSize = 0.5) + scale_y_continuous(limits = c(0, 13))
gwas_bwt_HAP_snp = ggman(gwas_emmax_snp, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt",
                       relative.positions = TRUE, title = "Hap EMMAx SNP", sigLine = -log10(2.6e-5),
                       pointSize = 0.5) + scale_y_continuous(limits = c(0, 13))
gwas_bwt_HAP_loco = ggman(gwas_emmax_loco, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt",
                        relative.positions = TRUE, title = "Hap EMMAx SNP LOCO", sigLine = -log10(2.6e-5),
                        pointSize = 0.5) + scale_y_continuous(limits = c(0, 13))
gwas_comp = plot_grid(gwas_bwt_GEMMA_ped, gwas_bwt_GEMMA_snp, gwas_bwt_GEMMA_loco,
          gwas_bwt_HAP_ped, gwas_bwt_HAP_snp, gwas_bwt_HAP_loco, ncol = 3)
save_plot("~/Dropbox/labbio/data/Atchley project/haplotype_inference_figures/gwas_comparison.png", gwas_comp, base_height = 5, base_aspect_ratio = 1.3, ncol = 3, nrow = 2)

founders_dict[[(gwas_rsnp_loco %>% filter(chr == 6) %>% arrange(p_lrt))[1,"rs"]]]

i = (gwas_emmax_loco %>% filter(chr == 6) %>% arrange(p_lrt))[1,"rs"]
i = which(gen$ID == i)
X = f6_model_matrices[,,i]
X_t <- crossprod(M[chr,,], X)
summary(lsfit(cbind(Xo_t[chr,,],X_t), Y_t[,chr], intercept = FALSE))

Y = f6_snped$Final_weight
Y[is.na(Y)] = mean(Y, na.rm = TRUE)
library(stanAnimal)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)
g = lmm_animal(Y, Xo, Af6_snp, chains = 4, iter = 2000, warmup = 1000)
g_ped = lmm_animal(Y, Xo, Af6, chains = 4, iter = 2000, warmup = 1000)
print(g, digits = 5)
V = sapply(rstan::extract(g, pars = c("sigma_G", "sigma_E")), mean)
V[1]/(V[1] + V[2])

library(pedigreemm)
ped17 <- pedigree(ped2$sire, ped2$dam, ped2$ID)  #restructure ped file
data = data.frame(Y = Y, sex = f6_snped$Sex, ID = f6_snped$ID)
mod_animalREML<-pedigreemm(Y ~ sex + (1|ID), pedigree=list(ID=ped17),
                           data = data, REML=TRUE,
                           control = lmerControl(check.nobs.vs.nlev="ignore",
                                                 check.nobs.vs.nRE="ignore"))
summary(mod_animalREML)
