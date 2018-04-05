source("./read_phenotypes.R")

raw_gen = read_table2("./data/filtered_atchley_genotypes_raw.txt", comment = "#")[-1]
tail(colnames(raw_gen))
n = length(colnames(raw_gen))

gen = dplyr::select(raw_gen, Chr_id, Start, everything()) %>%
  rename(chr = Chr_id,
         pos = Start) %>%
  arrange(chr, pos) %>%
  rename(ID = Affy_SNP_ID) %>%
  filter(chr != 21)

## From http://churchill-lab.jax.org/mousemapconverter
## NCBI Build 37 bp -> Sex-Averaged cM Cox
positions = read.table("./data/marker_positions.txt", header = FALSE, stringsAsFactors = FALSE) %>% tbl_df
names(positions) = c("chr2", "pos", "chr", "gpos")
positions$chr = as.integer(positions$chr)
tail(positions$chr)
gen = inner_join(positions, gen, by = c("chr", "pos")) %>%
  dplyr::select(ID, chr, everything(), -chr2)

IDs = data_frame(ID = strsplit(colnames(raw_gen)[c(-(n-2), -(n-1), -n)], split = "_") %>%
                   map(1) %>%
                   unlist %>%
                   as.numeric,
                 pID = colnames(raw_gen)[c(-(n-2), -(n-1), -n)])

full_snped = inner_join(full_data, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID, Final_weight)

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)

f5_snped = inner_join(full_data_F5, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)

f1_snped = inner_join(full_data_F1, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)

pat_snped = inner_join(full_data_Strain, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)

nrow(pat_snped)

