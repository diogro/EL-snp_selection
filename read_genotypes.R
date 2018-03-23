library(readr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(doMC)
n_chunks = 8
registerDoMC(n_chunks)

line_order = c("A13", "A31", "A41", "A23", "A22", "A42")[6:1]

full_data_F6 = read_csv("./data/Mouse phenotypes.csv") %>%
  dplyr::select(Litter_ID_new:Sex,
                Gen, Pat_ID, Mat_ID, Nurse_ID, Litter_size_birth,
                Birth_litter_size_weaning, Foster_litter_size_weaning,
                Weight_D0:Weight_D70, Final_weight, Liver:Fat) %>%
  filter(Gen == "F6")

full_data_F5 = read_csv("./data/Mouse phenotypes.csv") %>%
  dplyr::select(Litter_ID_new:Sex,
                Gen, Pat_ID, Mat_ID, Nurse_ID, Litter_size_birth,
                Birth_litter_size_weaning, Foster_litter_size_weaning,
                Weight_D0:Weight_D70, Final_weight, Liver:Fat) %>%
  filter(Gen == "F5")

full_data_F1 = read_csv("./data/Mouse phenotypes.csv") %>%
  dplyr::select(Litter_ID_new:Sex,
                Gen, Pat_ID, Mat_ID, Nurse_ID, Litter_size_birth,
                Birth_litter_size_weaning, Foster_litter_size_weaning,
                Weight_D0:Weight_D70, Final_weight, Liver:Fat) %>%
  filter(Gen == "F1")

full_data_Strain = read_csv("./data/Mouse phenotypes.csv") %>%
  dplyr::select(Litter_ID_new:Sex,
                Gen, Strain, Pat_ID, Mat_ID, Nurse_ID, Litter_size_birth,
                Birth_litter_size_weaning, Foster_litter_size_weaning,
                Weight_D0:Weight_D70, Final_weight, Liver:Fat) %>%
  filter(Gen == "Strain")

full_data_F6$ID[full_data_F6$ID == 3202] = 3302

pedigree = as.data.frame(read.csv("./data/Intercross_pedigree.csv")) %>%
  rename(id = animal) %>% orderPed

raw_gen = read_table2("~/Dropbox/labbio/data/Atchley project/Genotypes/raw_atchley_genotypes.txt", comment = "#")[-1]
tail(colnames(raw_gen))
n = length(colnames(raw_gen))

gen = dplyr::select(raw_gen, Chr_id, Start, everything()) %>%
  select(-Allele_Count) %>%
  rename(chr = Chr_id,
         pos = Start) %>%
  arrange(chr, pos) %>%
  filter(chr != 21)

## From http://churchill-lab.jax.org/mousemapconverter
## NCBI Build 37 bp -> Sex-Averaged cM Cox
positions = read_table2("./data/marker_positions.txt", col_names = c("chr", "pos", "chr2", "gpos"), col_types = "iidd")[,-3]

gen = inner_join(positions, gen, by = c("chr", "pos"))

IDs = data_frame(ID = strsplit(colnames(raw_gen)[c(-(n-2), -(n-1), -n)], split = "_") %>% map(1) %>% unlist %>% as.numeric,           pID = colnames(raw_gen)[c(-(n-2), -(n-1), -n)])

pat_snped = inner_join(full_data_Strain, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)
table(pat_snped$Strain)

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)

f5_snped = inner_join(full_data_F5, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)

f1_snped = inner_join(full_data_F1, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)

anti_join(IDs, pat_snped, by = "ID") %>%
  anti_join(f6_snped, by = "ID") %>%
  anti_join(f5_snped, by = "ID") %>%
  anti_join(f1_snped, by =  "ID")

nrow(IDs) - (nrow(pat_snped) + nrow(f6_snped) + nrow(f5_snped) + nrow(f1_snped))

pat_gen = gen %>% select(chr, pos, gpos, pat_snped$pID)
x = pat_gen[7,-c(1, 2, 3)]
(t1 = table(pat_snped$Strain[x == "AA"]))
(t2 = table(pat_snped$Strain[x == "BB"]))
(t3 = table(pat_snped$Strain[x == "AB"]))

toNumeric  = function(x) factor(x, levels = c('AA', 'AB', 'BB'))

current_line = line_order[3]
current_line_gen = gen %>% select(Chr_id, Start, filter(pat_snped, Strain == current_line)$pID)

A23 = map_df(current_line_gen[,-c(1, 2)], toNumeric) %>% as.data.frame %>% na.omit
head(A23)
cov(A23)

gen_class = c("AA", "BB", "AB")
cl = line_order[1]
gen_class[c(cl %in% names(t1), cl %in% names(t2), cl %in% names(t3))]


is_any_kind = function(x){
  tx = table(x)
  snp_kind = NULL
  if(!is.null(snp_kind <- get_private(x, tx))) return(c("private", snp_kind))
  if(!is.null(snp_kind <- get_two_four(x, tx))) return(c("two_four", snp_kind))
  if(!is.null(snp_kind <- get_50_50(x, tx))) return(c("50_50", snp_kind))
  if(any(tx == 1)){
    if(get_A22_p_het(x)) return(c("private_het", "A22_pri_het"))
  } else if(any(tx == 3)){
    if(get_A22_het_p_42(x)) return(c("private_het", "A22_het_p_A42"))
    else if(get_A22_het_p_23(x)) return(c("private_het", "A22_het_p_A23"))
  }
  return(c(NA, NA))
}
