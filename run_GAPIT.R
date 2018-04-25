#Step 0: Import library and GAPIT functions run this section each time to start R)
#######################################################################################
library('MASS') # required for ginv
if(!require(multtest)){
  source("http://www.bioconductor.org/biocLite.R");
  biocLite("multtest");
  library(multtest)}
if(!require(gplots)){install.packages("gplots"); library(gplots)}
if(!require(compiler)){install.packages("compiler"); library(compiler)}
if(!require(scatterplot3d)){install.packages("scatterplot3d"); library(scatterplot3d)}
if(!require(data.table)){install.packages("data.table"); library(data.table)}

source("http://www.zzlab.net/GAPIT/emma.txt")
source("/home/diogro/projects/EL-snp_selection/gapit_functions.R")

source("read_genotypes.R")

f5f6_snped = inner_join(full_data_F5F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID, contains("growth"), Final_weight)

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID, contains("growth"), Final_weight) %>% na.omit

phenotypes = select(f6_snped, ID, Final_weight, contains("growth")) %>% distinct(ID, .keep_all = TRUE) %>% rename(taxa = ID)
covars = select(f6_snped, ID, Sex) %>% distinct(ID, .keep_all = TRUE) %>% mutate(Sex = .n(Sex) - 1) %>% rename(taxa = ID)

f6_gen = gen %>% dplyr::select(ID, chr, pos, gpos, as.character(f6_snped$ID))
dt_gen = data.table(f6_gen[,-c(1:4)])
for(col in names(dt_gen)) set(dt_gen, i=which(dt_gen[[col]]=="1/1"), j=col, value=0)
for(col in names(dt_gen)) set(dt_gen, i=which(dt_gen[[col]]=="0/1"), j=col, value=1)
for(col in names(dt_gen)) set(dt_gen, i=which(dt_gen[[col]]=="0/0"), j=col, value=2)
dt_gen <- dt_gen[, lapply(.SD, as.numeric)]
qtl_rel_gen = t(dt_gen)
dimnames(qtl_rel_gen)[[2]] = f6_gen$ID
myGD = data.frame(taxa = f6_snped$ID, qtl_rel_gen)
myGD[1:10, 1:10]
map = f6_gen[,1:3]

setwd("/home/diogro/projects/EL-snp_selection/data/GAPIT/no_group")
myGAPIT_no_group <- GAPIT(
  Y=as.data.frame(phenotypes),
  CV=as.data.frame(covars),
  GD=myGD,
  GM=as.data.frame(map),
  PCA.total = 40,
  group.by = 1,
  group.from = 1274,
  group.to = 1274
)

setwd("/home/diogro/projects/EL-snp_selection/data/GAPIT/auto_grouped/")
myGAPIT_no_group <- GAPIT(
  Y=as.data.frame(phenotypes),
  CV=as.data.frame(covars),
  GD=myGD,
  GM=as.data.frame(map),
  PCA.total = 40
)

setwd("/home/diogro/projects/EL-snp_selection/data/GAPIT/SUPER/")
myGAPIT_no_group <- GAPIT(
  Y=as.data.frame(phenotypes),
  CV=as.data.frame(covars),
  GD=myGD,
  GM=as.data.frame(map),
  PCA.total = 40,
  LD = 0.01,
  sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER
  sangwich.bottom="SUPER" #options are GLM,MLM,CMLM, FaST and SUPER
)

