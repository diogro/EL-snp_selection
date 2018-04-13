library(data.table)

raw_gen = read_table2("./data/filtered_atchley_genotypes_raw.txt", comment = "#")[-1]
tail(colnames(raw_gen))
n = ncol(raw_gen)
IDs = data_frame(pID = (colnames(raw_gen)[-c(n-2, n-1, n)]))

IDs = data_frame(ID = strsplit(colnames(raw_gen)[c(-(n-2), -(n-1), -n)], split = "_") %>%
                 map(1) %>%
                 unlist %>%
                 as.numeric,
                 pID = colnames(raw_gen)[c(-(n-2), -(n-1), -n)])

full_snped = inner_join(full_data,
                        IDs,
                        by = "ID") %>%
  select(ID, pID, Litter_ID_new:Mat_ID, Final_weight)

plink_map_raw = read_table2("./data/plink_files/filtered_atchley_genotypes_plink_raw.map",
                            comment = "#", col_names = c("chr", "ID", "gpos", "pos"))

plink_map =
  inner_join(select(plink_map_raw, -gpos),
             positions,
             by = c("chr", "pos")) %>%
  mutate(gpos = gpos / 100) %>%
  dplyr::select(chr, ID, gpos, pos)

write_delim(plink_map, "./data/plink_files/filtered_atchley_genotypes_plink.map", col_names = FALSE, delim = " ")

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
  dplyr::select("Family", "ID", "father", "mother", "sex", "phenotype", everything(), -pID) %>%
  distinct(ID, .keep_all = TRUE)

write_delim(new_ped_file, "./data/plink_files/filtered_atchley_genotypes_plink.ped", col_names = FALSE, delim = " ")
system("sed -i 's/\"//g' ./data/plink_files/filtered_atchley_genotypes_plink.ped")

new_pheno_file = select(new_ped_file, Family, ID, phenotype) %>% na.omit
write_delim(new_pheno_file, "./data/fastlmm/fastlmm_pheno.txt", col_names = FALSE,delim = " ")

new_cov_file = select(new_ped_file, Family, ID, sex) %>% mutate(sex = sex -1)
write_delim(new_cov_file, "./data/fastlmm/fastlmm_cov.txt", col_names = FALSE,delim = " ")
