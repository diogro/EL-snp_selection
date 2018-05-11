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
if(!require(qvalue)){source("https://bioconductor.org/biocLite.R"); biocLite("qvalue"); library(qvalue)}
if(!require(ggman)){devtools::install_github("drveera/ggman"); library(ggman)}

n_chunks = 6
registerDoMC(n_chunks)

line_order = c("A13", "A31", "A41", "A23", "A22", "A42")[6:1]

full_data = read_csv("./data/Mouse phenotypes.csv") %>%
    mutate(ID = as.character(ID)) %>%
    dplyr::select(ID, Litter_ID_new:Sex,
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
full_data$ID = as.numeric(full_data$ID)
ped2 = (left_join(dplyr::rename(pedigree, ID = id), dplyr::select(full_data, ID, Sex), by = "ID"))
full_data$ID = as.character(full_data$ID)
missing = full_snped[!full_snped$ID %in% ped2$ID,c("ID", "Mat_ID", "Pat_ID", "Sex")]
names(missing) = names(ped2)
ped2 = rbind(ped2, missing) %>% orderPed
missing_sire = data.frame(ID = unique(na.omit(ped2$sire[!ped2$sire %in% ped2$ID])), dam = NA, sire = NA, Sex = "M")
ped2 = rbind(ped2, missing_sire) %>% orderPed
missing_dam = data.frame(ID = unique(na.omit(ped2$dam[!ped2$dam %in% ped2$ID])), dam = NA, sire = NA, Sex = "F")
ped2 = rbind(ped2, missing_dam) %>% orderPed
for(i in 1:nrow(ped2)){
  if(is.na(ped2$Sex[i])){
    if(ped2$ID[i] %in% ped2$dam) ped2$Sex[i] = "F"
    else if(ped2$ID[i] %in% ped2$sire) ped2$Sex[i] = "M"
    else ped2$Sex[i] = "M"
  }
}

A = 2*kinship(pedigree)
A = A + diag(nrow(A)) * 1e-4
ids = as.character(f6_snped$ID)
Af6 = (A[ids,ids])
colnames(Af6) = rownames(Af6) = f6_snped$ID



#source("./create_kinship_king.R")
