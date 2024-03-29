if(!require(lme4)){install.packages("lme4"); library(lme4)}
if(!require(lme4qtl)){devtools::install_github("variani/lme4qtl"); library(lme4qtl)}
if(!require(lmerTest)){devtools::install_github("runehaubo/lmerTest"); library(lmerTest)}
if(!require(stanAnimal)){devtools::install_github("diogro/stanAnimal", subdir = "package"); library(stanAnimal)}
if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(doMC)){install.packages("doMC"); library(doMC)}
if(!require(MCMCglmm)){install.packages("MCMCglmm"); library(MCMCglmm)}
# if(!require(MasterBayes)){install.packages("MasterBayes"); library(MasterBayes)}
if(!require(ggplot2)){devtools::install_github("tidyverse/ggplot2"); library(ggplot2)}
if(!require(cowplot)){devtools::install_github("wilkelab/cowplot"); library(cowplot)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}
if(!require(viridis)){install.packages("viridis"); library(viridis)}
if(!require(QTLRel)){install.packages("QTLRel"); library(QTLRel)}
if(!require(qvalue)){source("https://bioconductor.org/biocLite.R"); biocLite("qvalue"); library(qvalue)}
if(!require(ggman)){devtools::install_github("drveera/ggman"); library(ggman)}
library(AtchleyMice)

registerDoMC(6)

line_order = c("A13", "A31", "A41", "A23", "A22", "A42")[6:1]

data(mice_info)

full_data = mice_info$full

#full_data <- mutate(full_data,
                    #growth_D0D7   = Weight_D7 - Weight_D0,
                    #growth_D7D14  = Weight_D14 - Weight_D7,
                    #growth_D14D21 = Weight_D21 - Weight_D14,
                    #growth_D21D28 = Weight_D28 - Weight_D21,
                    #growth_D28D35 = Weight_D35 - Weight_D28,
                    #growth_D35D42 = Weight_D42 - Weight_D35,
                    #growth_D42D49 = Weight_D49 - Weight_D42,
                    #growth_D49D56 = Weight_D56 - Weight_D49)

# full_data <- mutate(full_data,
#                     growth_D0D14  = Weight_D14 - Weight_D0,
#                     growth_D14D28 = Weight_D28 - Weight_D14,
#                     growth_D28D42 = Weight_D42 - Weight_D28,
#                     growth_D42D56 = Weight_D56 - Weight_D42)

full_data_F5F6 =  full_data %>% filter(Gen == "F6" | Gen == "F5")

full_data_F6 =  full_data %>% filter(Gen == "F6")

full_data_F5 = full_data %>% filter(Gen == "F5")

full_data_F1 = full_data %>% filter(Gen == "F1")

full_data_Strain = full_data %>% filter(Gen == "Strain")

# growth_traits = c("growth_D0D14", "growth_D14D28", "growth_D28D42", "growth_D42D56")

# full_data_F6[,growth_traits]$growth_D42D56
#source("./create_kinship_king.R")
