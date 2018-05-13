source("./read_phenotypes.R")

#raw_gen = read_table2("./data/filtered_atchley_genotypes_raw.txt", comment = "#")[-1]
#raw_gen2 = read_table2("./data/plink_files/missing.vcf", comment = "##") %>% arrange(`#CHROM`, POS)
raw_gen = read_table2("./data/plink_files/atchley_imputed_text.vcf", comment = "##") %>% arrange(`#CHROM`, POS)

colnames(raw_gen)[-c(1:9)] = as.character(as.numeric(colnames(raw_gen)[-c(1:9)]))

tail(colnames(raw_gen))
n = length(colnames(raw_gen))

gen = dplyr::select(raw_gen, `#CHROM`, POS, ID, everything()) %>%
  dplyr::select(-REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT) %>%
  rename(chr = `#CHROM`,
         pos = POS) %>%
  arrange(chr, pos) %>%
  filter(chr != 21)

## From http://churchill-lab.jax.org/mousemapconverter
## NCBI Build 37 bp -> Sex-Averaged cM Cox
positions = read.table("./data/marker_positions.txt", header = FALSE, stringsAsFactors = FALSE) %>% tbl_df
names(positions) = c("chr2", "pos", "chr", "gpos")
positions$chr = as.integer(positions$chr)
tail(positions$chr)
gen = inner_join(positions, gen, by = c("chr", "pos")) %>%
  dplyr::select(ID, chr, everything(), -chr2)

IDs = data_frame(ID = as.character(colnames(raw_gen)[-c(1:9)]))

full_snped = inner_join(full_data, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID, Final_weight)

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID, Final_weight)

f5_snped = inner_join(full_data_F5, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID)

f1_snped = inner_join(full_data_F1, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID)

pat_snped = inner_join(full_data_Strain, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID)

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
load("./data/gemma_relatedness.Rdata")
#pedAll <- pedigree(id=ped2$ID,
                   #dadid=ped2$sire, momid=ped2$dam,
                   #sex=ped2$Sex)

