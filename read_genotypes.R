library(readr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)
library(tidyr)
if(!require(purrr)){install.packages("purrr"); library(purrr)}
library(doMC)
library(MasterBayes)
library(dplyr)
n_chunks = 4
registerDoMC(n_chunks)

line_order = c("A13", "A31", "A41", "A23", "A22", "A42")[6:1]

full_data = read_csv("./data/Mouse phenotypes.csv") %>%
  dplyr::select(Litter_ID_new:Sex,
                Gen, Pat_ID, Mat_ID, Nurse_ID, Strain, Litter_size_birth,
                Birth_litter_size_weaning, Foster_litter_size_weaning,
                Weight_D0:Weight_D70, Final_weight, Liver:Fat)

full_data_F6 =  full_data %>% filter(Gen == "F6")

full_data_F5 = full_data %>% filter(Gen == "F5")

full_data_F1 = full_data %>% filter(Gen == "F1")

full_data_Strain = full_data %>% filter(Gen == "Strain")

full_data_F6$ID[full_data_F6$ID == 3202] = 3302

pedigree = as.data.frame(read.csv("./data/Intercross_pedigree.csv")) %>%
  rename(id = animal) %>% orderPed

raw_gen = read_table2("./data/raw_atchley_genotypes.txt", comment = "#")[-1]
tail(colnames(raw_gen))
n = length(colnames(raw_gen))

gen = dplyr::select(raw_gen, Chr_id, Start, everything()) %>%
  dplyr::select(-Allele_Count) %>%
  rename(chr = Chr_id,
         pos = Start) %>%
  arrange(chr, pos) %>%
  mutate(ID = paste(chr, pos, sep = "_")) %>%
  filter(chr != 21)

## From http://churchill-lab.jax.org/mousemapconverter
## NCBI Build 37 bp -> Sex-Averaged cM Cox
positions = read.table("./data/marker_positions.txt", header = FALSE, stringsAsFactors = FALSE) %>% tbl_df
names(positions) = c("chr2", "pos", "chr", "gpos")
positions$chr = as.integer(positions$chr)
tail(positions$chr)
gen = inner_join(positions, gen, by = c("chr", "pos")) %>% dplyr::select(ID, chr, everything(), -chr2)

IDs = data_frame(ID = strsplit(colnames(raw_gen)[c(-(n-2), -(n-1), -n)], split = "_") %>%
                   map(1) %>%
                   unlist %>%
                   as.numeric,
                 pID = colnames(raw_gen)[c(-(n-2), -(n-1), -n)])

pat_snped = inner_join(full_data_Strain, IDs, by = "ID") %>%
  dplyr::select(pID, ID, Strain, Litter_ID_new:Mat_ID)
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

pat_gen = gen %>% dplyr::select(ID, chr, pos, gpos, pat_snped$pID)
ggplot(pat_gen, aes(x = gpos, y = chr)) + geom_point(size = 0.3, alpha = 0.5)
current_line = line_order[2]
getConsensusCall = function(current_line){
  current_line_gen = gen %>% dplyr::select(ID, chr, pos, gpos, filter(pat_snped, Strain == current_line)$pID)
  getConsensusCallPerSnp = function(current_snp){
    x = current_snp[-c(1, 2, 3, 4)]
    tx = table(x)
    colnames(tx)[which.max(tx)]
  }
  unlist(alply(current_line_gen, 1, getConsensusCallPerSnp))
}
strain_genotypes = laply(line_order, getConsensusCall, .progress = "text")
rownames(strain_genotypes) = line_order

header = paste("marker", dim(strain_genotypes)[2], "strain 6")
strain_string = paste0("strain_names ", paste(line_order, collapse = " "))
NA_string = paste0("allele NA ", paste(rep(1/6, 6), collapse = " "))
generateSNPentry = function(i){
  current_snp_ID = pat_gen[i,1:4]
  current_snp = strain_genotypes[,i]
  (hasAA = current_snp == "AA")
  (hasBB = current_snp == "BB")
  (hasAB = current_snp == "AB")
  hasA = (map2_lgl(hasAA, hasAB, `|`))
  hasB = (map2_lgl(hasBB, hasAB, `|`))
  nA = sum(hasA)
  nB = sum(hasB)
  rowA = numeric(6)
  rowB = numeric(6)
  rowA[hasA] = 1/nA
  rowB[hasB] = 1/nB
  out = paste0("marker ", current_snp_ID$ID, " 3 ", current_snp_ID$gpos, "\n",
               NA_string, "\n",
               "allele A ", paste(rowA, collapse = " "), "\n",
               "allele B ", paste(rowB, collapse = " "))
  return(out)
}
marker_entries = llply(seq_along(1:dim(strain_genotypes)[2]), generateSNPentry, .progress = "text")
marker_entries[10000]

write_lines(header, "./data/sink_alleles.txt")
write_lines(strain_string, "./data/sink_alleles.txt", append = TRUE)
write_lines(marker_entries, "./data/sink_alleles.txt", append = TRUE)


f6_genotypes = gen %>% select(ID, chr, pos, gpos, filter(f6_snped)$pID)

df = f6_genotypes[,-c(2:4)] %>%
  replace(., . == "NoCall", "NA\tNA") %>%
  replace(., . == "AA", "A\tA") %>%
  replace(., . == "BB", "B\tB") %>%
  replace(., . == "AB", "A\tB")

t.df = df %>%
  gather(var, value, -ID) %>%
  spread(ID, value)
t.df = t.df %>%
  rename(pID = var)

f6_happy = inner_join(dplyr::select(inner_join(dplyr::select(full_data_F6, ID, Final_weight), f6_snped, by = "ID"), pID, Final_weight), t.df, by = "pID") %>%
  rename(SAMPLE_ID = pID, PHENOTYPE = Final_weight) %>%
  select(SAMPLE_ID, PHENOTYPE, df$ID)

f6_happy

write_tsv(f6_happy, "./data/f6_genotypes_happy.tsv")


