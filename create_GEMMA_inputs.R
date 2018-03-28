source("read_genotypes.R")


gen6 = gen %>% filter(chr == 6)

f5f6_snped = inner_join(full_data_F5F6, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)

f5_snped = inner_join(full_data_F5, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)

f6_genotypes = gen6 %>% select(ID, chr, pos, gpos, filter(f6_snped)$pID)

df = f6_genotypes %>%
  replace(., . == "NoCall", NA) %>%
  replace(., . == "AA", "0") %>%
  replace(., . == "BB", "1") %>%
  replace(., . == "AB", "0.5")

for(i in 1:dim(df)[1]){
  missing = which(is.na(df[i,]))
  if(length(missing) > 0)
    for(j in 1:length(missing)){
      ind = colnames(df)[missing[j]]
      litter = f6_snped$Litter_ID_new[f6_snped$pID == ind]
      df[i, missing] = mean(as.numeric(df[missing[j],filter(f6_snped, Litter_ID_new == litter)$pID]), na.rm = TRUE)
    }
}

