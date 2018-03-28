source("./read_genotypes.R")

gen = gen %>% filter(chr == 6)

IDs = data_frame(ID = strsplit(colnames(raw_gen)[c(-(n-2), -(n-1), -n)], split = "_") %>%
                   map(1) %>%
                   unlist %>%
                   as.numeric,
                 pID = colnames(raw_gen)[c(-(n-2), -(n-1), -n)])

pat_snped = inner_join(full_data_Strain, IDs, by = "ID") %>%
  dplyr::select(pID, ID, Strain, Litter_ID_new:Mat_ID)

f5f6_snped = inner_join(full_data_F5F6, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)

f5_snped = inner_join(full_data_F5, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)

f1_snped = inner_join(full_data_F1, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)

pat_gen = gen %>%
  dplyr::select(ID, chr, pos, gpos, pat_snped$pID)
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

header = paste("markers", dim(strain_genotypes)[2], "strains 6")
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
marker_entries[100]

write_lines(header, "./data/markers_strain_split_chr6.txt")
write_lines(strain_string, "./data/markers_strain_split_chr6.txt", append = TRUE)
write_lines(marker_entries, "./data/markers_strain_split_chr6.txt", append = TRUE)

f5f6_genotypes = gen %>% select(ID, chr, pos, gpos, filter(f6_snped)$pID, filter(f5_snped)$pID)

df = f5f6_genotypes[,-c(2:4)] %>%
  replace(., . == "NoCall", "NA\tNA") %>%
  replace(., . == "AA", "A\tA") %>%
  replace(., . == "BB", "B\tB") %>%
  replace(., . == "AB", "A\tB")

t.df = df %>%
  gather(var, value, -ID) %>%
  spread(ID, value)
t.df = t.df %>%
  rename(pID = var)

f5f6_happy =
  inner_join(dplyr::select(inner_join(dplyr::select(full_data_F5F6, ID, Final_weight), f5f6_snped, by = "ID"),                                    Litter_ID_new, ID, Mat_ID, Pat_ID, Sex, pID, Final_weight),
             t.df, by = "pID") %>%
  rename("#Family-id" = Litter_ID_new,
         "individual-id" = ID,
         "mother-id" = Mat_ID,
         "father-id" = Pat_ID,
         sex = Sex,
         phenotype = Final_weight) %>%
  mutate(sex = as.numeric(as.factor(sex))) %>%
  dplyr::select("#Family-id", "individual-id", "mother-id", "father-id", "sex", "phenotype", df$ID)

f6_happy

write_tsv(f5f6_happy, "./data/f5f6_genotypes_happy_chr6.PED")

system("sed -i 's/\"//g' ./data/f5f6_genotypes_happy_chr6.PED")