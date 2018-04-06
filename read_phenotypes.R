if(!require(lme4)){install.packages("lme4"); library(lme4)}
if(!require(lme4qtl)){devtools::install_github("variani/lme4qtl"); library(lme4qtl)}
if(!require(lmerTest)){devtools::install_github("runehaubo/lmerTest"); library(lmerTest)}
if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(doMC)){install.packages("doMC"); library(doMC)}
if(!require(MCMCglmm)){install.packages("MCMCglmm"); library(MCMCglmm)}
if(!require(MasterBayes)){install.packages("MasterBayes"); library(MasterBayes)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}
if(!require(viridis)){install.packages("viridis"); library(viridis)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}
if(!require(QTLRel)){install.packages("QTLRel"); library(QTLRel)}

n_chunks = 6
registerDoMC(n_chunks)

line_order = c("A13", "A31", "A41", "A23", "A22", "A42")[6:1]

full_data = read_csv("./data/Mouse phenotypes.csv") %>%
  dplyr::select(Litter_ID_new:Sex,
                Gen, Pat_ID, Mat_ID, Nurse_ID, Strain, Litter_size_birth,
                Birth_litter_size_weaning, Foster_litter_size_weaning,
                Weight_D0:Weight_D70, Final_weight, Liver:Fat)
full_data$ID[full_data$ID == 3202] = 33020

full_data <- mutate(full_data,
                    growth_D0D3   = Weight_D3 - Weight_D0,
                    growth_D3D7   = Weight_D7 - Weight_D3,
                    growth_D7D14  = Weight_D14 - Weight_D7,
                    growth_D14D21 = Weight_D21 - Weight_D14,
                    growth_D21D28 = Weight_D28 - Weight_D21,
                    growth_D28D35 = Weight_D35 - Weight_D28,
                    growth_D35D42 = Weight_D42 - Weight_D35,
                    growth_D42D49 = Weight_D49 - Weight_D42,
                    growth_D49D56 = Weight_D56 - Weight_D49)

full_data_F5F6 =  full_data %>% filter(Gen == "F6" | Gen == "F5")

full_data_F6 =  full_data %>% filter(Gen == "F6")

full_data_F5 = full_data %>% filter(Gen == "F5")

full_data_F1 = full_data %>% filter(Gen == "F1")

full_data_Strain = full_data %>% filter(Gen == "Strain")

.n = function(x) as.numeric(factor(x, levels = c("M", "F")))
pedigree = as.data.frame(read.csv("./data/Intercross_pedigree2.csv")) %>% dplyr::rename(id = animal) %>% orderPed
ped2 = (left_join(dplyr::rename(pedigree, ID = id), dplyr::select(full_data, ID, Sex), by = "ID"))
for(i in 1:nrow(ped2)){
  if(is.na(ped2$Sex[i])){
    if(ped2$ID[i] %in% ped2$dam) ped2$Sex[i] = "F"
    else if(ped2$ID[i] %in% ped2$sire) ped2$Sex[i] = "M"
    else ped2$Sex[i] = "M"
  }
}
pedAll <- pedigree(id=ped2$ID,
                   dadid=ped2$sire, momid=ped2$dam,
                   sex=ped2$Sex)

