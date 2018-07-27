source("read_phenotypes.R")
source("read_genotypes.R")
source("data/LASSO_cluster/emma.r")

if(!require(foreach)){install.packages("foreach"); library(foreach)}
if(!require(doMC)){install.packages("doMC"); library(doMC)}
registerDoMC(2)

load("./data/f6_model_matrices.Rdata")

f6_snped = f6_snped[match(dimnames(f6_model_matrices)[[1]], f6_snped$ID),]
all(f6_snped$ID == dimnames(f6_model_matrices)[[1]])

Xo <- cbind(1, .n(f6_snped$Sex)-1,
            scale(cbind(f6_snped$Litter_size_birth,
                  f6_snped$Birth_litter_size_weaning, f6_snped$Foster_litter_size_weaning)))
colnames(Xo) = c("Intercept", "sex", "lsb", "blsw", "flsw")


X = f6_model_matrices[,,1]
lm_no_int = lm(scale(Final_weight)~ 0 + cbind(Xo,X), data = f6_snped)
lm_int = lmer(scale(Final_weight)~  0 + cbind(Xo,X) + (1|Mat_ID), data = f6_snped)
c_no_int = summary(lm_no_int)
c_no_int$coefficients[5:10,1]
c_int = summary(lm_int)
c_int$coefficients[6:10,1]

Af6 = A[f6_snped$ID, f6_snped$ID]
Af6_snp = A_snp[f6_snped$ID, f6_snped$ID]
Af6_loco = A_loco[,f6_snped$ID, f6_snped$ID]

source("hapGWAS.R")
Y = f6_snped$Final_weight
K_norm = Af6

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
gwas_comp = plot_grid(gwas_bwt_HAP_ped, gwas_bwt_HAP_snp, gwas_bwt_HAP_loco, ncol = 3)
save_plot("~/Dropbox/labbio/data/Atchley project/haplotype_inference_figures/gwas_comparison.png", gwas_comp, base_height = 5, base_aspect_ratio = 1.3, ncol = 3, nrow = 2)

founders_dict[[(gwas_rsnp_loco %>% filter(chr == 6) %>% arrange(p_lrt))[1,"rs"]]]

i = (gwas_emmax_loco %>% filter(chr == 6) %>% arrange(p_lrt))[1,"rs"]
i = which(gen$ID == i)
X = f6_model_matrices[,,i]
X_t <- crossprod(M[chr,,], X)
summary(lsfit(cbind(Xo_t[chr,,],X_t), Y_t[,chr], intercept = FALSE))

Y = f6_snped$Final_weight
Y[is.na(Y)] = mean(Y, na.rm = TRUE)
options(mc.cores = 4)
g = lmm_animal(Y, cbind(1, Xo[,1]), Af6_snp)
g_ped = lmm_animal(Y, cbind(1, Xo), Af6)
print(g, digits = 5)
V = sapply(rstan::extract(g, pars = c("sigma_G", "sigma_E")), mean)
V[1]/(V[1] + V[2])

install.packages("pedigreemm")
library(pedigreemm)
ped17 <- pedigree(ped2$sire, ped2$dam, ped2$ID)  #restructure ped file
data = data.frame(Y = Y, sex = f6_snped$Sex, ID = f6_snped$ID)
mod_animalREML<-pedigreemm(Y ~ sex + (1|ID), pedigree=list(ID=ped17),
                           data = data, REML=TRUE,
                           control = lmerControl(check.nobs.vs.nlev="ignore",
                                                 check.nobs.vs.nRE="ignore"))
summary(mod_animalREML)
