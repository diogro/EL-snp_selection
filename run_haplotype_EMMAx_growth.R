source("read_phenotypes.R")
source("read_genotypes.R")
source("data/LASSO_cluster/emma.r")

load("./data/f6_model_matrices.Rdata")

Xo <- cbind(rep(1, nrow(f6_snped)),.n(f6_snped$Sex)-1)
n = nrow(Af6)
#K_norm<-(n-1)/sum((diag(n)-matrix(1,n,n)/n)*Af6)*Af6
Af6_snp = A_snp[f6_snped$ID, f6_snped$ID]
Af6_loco = A_loco[,f6_snped$ID, f6_snped$ID]

hapGWAS <- function(Y, K_norm) {
  Y[is.na(Y)] = mean(Y, na.rm = TRUE)
  null<-emma.REMLE(Y,as.matrix(Xo),K_norm)

  #EMMAX

  M<-solve(chol(null$vg*K_norm+null$ve*diag(n)))
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
hist(pval, breaks = 200)
plot.inflation(pval)
plot(-log10(pval), pch = 19)

gwas_emmax_ped = hapGWAS(f6_snped$Final_weight, Af6)
gwas_emmax_snp = hapGWAS(f6_snped$Final_weight, Af6_snp)
gwas_emmax_loco = hapGWASLOCO(f6_snped$Final_weight, Af6_loco)
save(gwas_emmax_snp, gwas_emmax_ped, gwas_emmax_loco, file = "./data/haplotype_gwas_bwt.Rdata")
gwas_bwt_p_ped = ggman(gwas_emmax_ped, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt",
                        relative.positions = TRUE, title = "Hap EMMAx PED", sigLine = -log10(2.6e-5),
                        pointSize = 1) + scale_y_continuous(limits = c(0, 13))
gwas_bwt_p_snp = ggman(gwas_emmax_snp, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt",
                        relative.positions = TRUE, title = "Hap EMMAx SNP", sigLine = -log10(2.6e-5),
                        pointSize = 1) + scale_y_continuous(limits = c(0, 13))
gwas_bwt_p_loco = ggman(gwas_emmax_loco, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt",
                        relative.positions = TRUE, title = "Hap EMMAx SNP LOCO", sigLine = -log10(2.6e-5),
                        pointSize = 1) + scale_y_continuous(limits = c(0, 13))
plot_grid(gwas_bwt_p_ped, gwas_bwt_p_snp, gwas_bwt_p_loco)
