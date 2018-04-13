source("read_genotypes.R")

library(devtools)
install_github("HannahVMeyer/PhenotypeSimulator",build_vignettes = TRUE)
library(PhenotypeSimulator)

vignette("UsagePhenotypeSimulator", package = "PhenotypeSimulator")
vignette("Simulation-and-LinearModel", package = "PhenotypeSimulator")

Simulation-and-LinearModel

# Set parameters
genVar <- 0.4
noiseVar <- 1- genVar
shared <- 0.6
independent <- 1 - shared
kinship <- getKinship(N=100, X=genotypes_sd, verbose = FALSE)
phenotype <- runSimulation(10, P = 1, genotypefile = "./data/plink_files/atchely_imputed_thinned",
                           format = "plink", cNrSNP = 30, genVar = genVar, h2s = h2s,
                           phi = 0.6, delta = 0.3, distBetaGenetic = "unif", mBetaGenetic = 0.5,
                           sdBetaGenetic = 1, NrFixedEffects = 1,
                           NrConfounders = c(1),
                           pIndependentConfounders = c(1),
                           distConfounders = c("cat_norm"),
                           catConfounders = c(2), pcorr = 0.8,
                           verbose = TRUE)
phenotype$rawComponents$genotypes
#> Set seed: 219453
#> The total noise var