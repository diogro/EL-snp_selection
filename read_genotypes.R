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

IDs = data_frame(ID = as.numeric(colnames(raw_gen)[-c(1:9)]))

full_snped = inner_join(full_data, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID, Final_weight)

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID)

f5_snped = inner_join(full_data_F5, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID)

f1_snped = inner_join(full_data_F1, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID)

pat_snped = inner_join(full_data_Strain, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID)

nrow(pat_snped)

#pedAll <- pedigree(id=ped2$ID,
                   #dadid=ped2$sire, momid=ped2$dam,
                   #sex=ped2$Sex)

