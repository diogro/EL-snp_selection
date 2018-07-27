hapGWAS <- function(Y, K_norm) {
  if(any(is.na(Y))) Y[is.na(Y)] = mean(Y, na.rm = TRUE)
  null<-emma.REMLE(Y, as.matrix(Xo), K_norm)

  #EMMAX
  n = nrow(K_norm)
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

  n = dim(K)[3]

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
    summary(model)
    pval[i] = comp$`Pr(>F)`[2]
  }
  close(pb)
  data.frame(rs = gen$ID, ps = gen$pos, chr = gen$chr, p_lrt = pval)
}
