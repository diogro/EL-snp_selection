source("read_phenotypes.R")
source("data/LASSO_cluster/emma.r")

load("./data/f6_model_matrices.Rdata")

null = lmer(Final_weight ~ Sex + (1|Mat_ID), data = f6_snped,  REML = FALSE)
p_values = vector("numeric", dim(f6_model_matrices)[[3]])
k = 1
for(i in k:length(p_values)){
  marker = lmer(Final_weight ~ Sex + f6_model_matrices[,,i] + (1|Mat_ID), data = f6_snped, REML = FALSE)
  aov = anova(null, marker)
  p_values[i] = aov$`Pr(>Chisq)`[[2]]
  k = i
}
qv = qvalue(p_values, fdr.level = 0.01)
hist(qv$qvalue, breaks = 500)
table(qv$significant)
plot(-log10(qv$qvalue), pch = 19)
plot(-log10(qv$qvalue), pch = 19)
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
plot.inflation(qv$qvalue)


A = 2*kinship(pedigree)
ids = as.character(f6_snped$ID)
Af6 = (A[ids,ids])
Af6 = Af6 + diag(nrow(Af6)) * 1e-4
colnames(Af6) = rownames(Af6) = f6_snped$ID
bend_Af6 = nearPD(as.matrix(Af6))
colnames(bend_Af6$mat) = rownames(bend_Af6$mat) = phenotypes$ID

Xo <- cbind(rep(1, nrow(f6_snped)),.n(f6_snped$Sex)-1)
n = nrow(Af6)
#K_norm<-(n-1)/sum((diag(n)-matrix(1,n,n)/n)*Af6)*Af6
K_norm = Af6
#NULL MODEL

Y = f6_snped$Final_weight
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
hist(pval, breaks = 200)
plot.inflation(pval)
plot(-log10(pval), pch = 19)

gwas_emmax = data.frame(rs = gen$ID, ps = gen$pos, chr = gen$chr, p_lrt = pval)
(gwas_bwt_p_ped = ggman(gwas_emmax, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_lrt",
                        relative.positions = TRUE, title = "Hap EMMAx", sigLine = -log10(2.6e-5),
                        pointSize = 1))
