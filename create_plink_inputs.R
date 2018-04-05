source("./read_genotypes.R")

library(data.table)

plink_map_raw = read_table2("./data/plink_files/filtered_atchley_genotypes_plink_raw.map", comment = "#", col_names = FALSE)
colnames(plink_map_raw) = c("chr", "ID", "gpos", "pos")

plink_map = inner_join(select(plink_map_raw, -gpos), positions, by = c("chr", "pos")) %>%
  dplyr::select(chr, ID, gpos, pos)

write_tsv(plink_map, "./data/plink_files/filtered_atchley_genotypes_plink.map", col_names = FALSE)

ped_file = tbl_df(fread("./data/plink_files/filtered_atchley_genotypes_plink_raw.ped", header = FALSE, skip = 1))
colnames(ped_file)[1] = "pID"
ped_file = mutate(ped_file, pID = paste0(pID, "_call_code"))
new_ped_file = inner_join(ped_file,
                          select(full_snped, Litter_ID_new, ID, pID, Mat_ID, Pat_ID, Sex, Final_weight),
                          by ="pID") %>%
  rename("Family" = Litter_ID_new,
         "mother" = Mat_ID,
         "father" = Pat_ID,
         sex = Sex,
         phenotype = Final_weight) %>%
  mutate(sex = as.numeric(factor(sex, levels = c("M", "F")))) %>%
  dplyr::select("Family", "ID", "father", "mother", "sex", "phenotype", everything(), -pID)

tb = table(new_ped_file$ID)
duplicates = names(tb)[tb > 1]
new_ped_file_non_dup = filter(new_ped_file, !ID %in% duplicates)

write_delim(new_ped_file_non_dup, "./data/filtered_atchley_genotypes_plink.ped", col_names = FALSE,delim = " ")
system("sed -i 's/\"//g' ./data/filtered_atchley_genotypes_plink.ped")

new_pheno_file$phenotype[is.na(new_pheno_file$phenotype)]
new_pheno_file = select(new_ped_file_non_dup, Family, ID, phenotype) %>% na.omit
write_delim(new_pheno_file, "./data/fastlmm_pheno.txt", col_names = FALSE,delim = " ")

new_cov_file = select(new_ped_file_non_dup, Family, ID, sex) %>% mutate(sex = sex -1)
write_delim(new_cov_file, "./data/fastlmm_cov.txt", col_names = FALSE,delim = " ")
