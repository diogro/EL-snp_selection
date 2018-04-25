source("./read_genotypes.R")

pat_snped = inner_join(full_data_Strain, IDs, by = "ID") %>%
  dplyr::select(ID, Strain, Litter_ID_new:Mat_ID)

f5f6_snped = inner_join(full_data_F5F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID)

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID)

f5_snped = inner_join(full_data_F5, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID)

f1_snped = inner_join(full_data_F1, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID)

pat_gen = gen %>% dplyr::select(ID, chr, pos, gpos, as.character(pat_snped$ID))

current_line = line_order[2]
getConsensusCall = function(current_line){
  current_line_gen = gen %>% dplyr::select(ID, chr, pos, gpos, as.character(filter(pat_snped, Strain == current_line)$ID))
  getConsensusCallPerSnp = function(current_snp){
    x = current_snp[-c(1, 2, 3, 4)]
    tx = table(x)
    colnames(tx)[which.max(tx)]
  }
  unlist(alply(current_line_gen, 1, getConsensusCallPerSnp))
}
strain_genotypes = laply(line_order, getConsensusCall, .progress = "text")
rownames(strain_genotypes) = line_order

header = paste("markers", dim(strain_genotypes)[2], "strains 6")
strain_string = paste0("strain_names ", paste(line_order, collapse = " "))
NA_string = paste0("allele NA ", paste(rep(1/6, 6), collapse = " "))
generateSNPentry = function(i){
  current_snp_ID = pat_gen[i,1:4]
  current_snp = strain_genotypes[,i]
  (hasAA = current_snp == "0/0")
  (hasBB = current_snp == "1/1")
  (hasAB = current_snp == "0/1")
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
               "allele 0 ", paste(rowA, collapse = " "), "\n",
               "allele 1 ", paste(rowB, collapse = " "))
  return(out)
}
marker_entries = llply(seq_along(1:dim(strain_genotypes)[2]), generateSNPentry, .progress = "text")
marker_entries[100]

write_lines(header, "./data/happy_markers_strain.txt")
write_lines(strain_string, "./data/happy_markers_strain.txt", append = TRUE)
write_lines(marker_entries, "./data/happy_markers_strain.txt", append = TRUE)

f6_genotypes = gen %>% select(ID, chr, pos, gpos, as.character(f6_snped$ID))

df = f6_genotypes[,-c(2:4)] %>%
  replace(., . == "./.", "NA\tNA") %>%
  replace(., . == "0/0", "0\t0") %>%
  replace(., . == "1/1", "1\t1") %>%
  replace(., . == "0/1", "0\t1")

t.df = df %>%
  gather(var, value, -ID) %>%
  spread(ID, value)
t.df = t.df %>%
  rename(ID = var) %>%
  mutate(ID = as.numeric(ID))

f6_happy =
  inner_join(dplyr::select(inner_join(dplyr::select(full_data_F5F6, ID, Final_weight), f5f6_snped, by = "ID"),                                    Litter_ID_new, ID, Mat_ID, Pat_ID, Sex, Final_weight),
             t.df, by = "ID") %>%
  rename("#Family-id" = Litter_ID_new,
         "individual-id" = ID,
         "mother-id" = Mat_ID,
         "father-id" = Pat_ID,
         sex = Sex,
         phenotype = Final_weight) %>%
  mutate(sex = as.numeric(as.factor(sex))) %>%
  dplyr::select("#Family-id", "individual-id", "mother-id", "father-id", "sex", "phenotype", df$ID)

write_tsv(f6_happy, "./data/happy_f6_genotypes.PED")
system("sed -i 's/\"//g' ./data/happy_f6_genotypes.PED")
