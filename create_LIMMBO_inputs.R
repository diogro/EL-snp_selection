source("read_genotypes.R")

f5f6_snped = inner_join(full_data_F5F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID)

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID, contains("growth"), Final_weight)

f5_snped = inner_join(full_data_F5, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID)

phenotypes = na.omit(select(f6_snped, ID, Sex, Final_weight, contains("growth")))

write_csv(select(phenotypes, ID, Final_weight, growth_D14D21), "./data/limmbo/limmbo_phenotypes.csv", col_names = TRUE)

#A = 2*snpKin
A = 2*kinship(pedigree)
ids = as.character(phenotypes$ID)
Af6 = (A[ids,ids])
colnames(Af6) = rownames(Af6) = phenotypes$ID
bend_Af6 = nearPD(as.matrix(Af6))
colnames(bend_Af6$mat) = rownames(bend_Af6$mat) = phenotypes$ID
write_csv(tbl_df(as.matrix(bend_Af6$mat)), "./data/limmbo/limmbo_relatedness.csv", col_names = TRUE)

write_csv(mutate(select(phenotypes, ID, Sex), Sex = .n(Sex)), "./data/limmbo/limmbo_covariates.csv", col_names = TRUE)
