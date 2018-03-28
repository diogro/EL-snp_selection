source("./read_genotypes.R")

ped_solar = ped2 %>%
  rename(id = ID,
         fa = sire,
         mo = dam,
         sex = Sex) %>%
  mutate(sex = .n(sex)) %>%
  select(id, fa, mo, sex)
ped_solar[is.na(ped_solar)] = 0
write_csv(ped_solar, "./data/solar_ped.ped")

gen = gen %>% filter(chr == 6)

f5f6_snped = inner_join(full_data_F5F6, IDs, by = "ID") %>%
  select(pID, ID, Litter_ID_new:Mat_ID)

f5f6_genotypes = gen %>% select(ID, chr, pos, gpos, filter(f6_snped)$pID, filter(f5_snped)$pID)

df = f5f6_genotypes[,-c(2:4)] %>%
  replace(., . == "NoCall", "-/-") %>%
  replace(., . == "AA", "A/A") %>%
  replace(., . == "BB", "B/B") %>%
  replace(., . == "AB", "A/B")

t.df = df %>%
  gather(var, value, -ID) %>%
  spread(ID, value)
t.df = t.df %>%
  rename(pID = var)

f5f6_solar =
  inner_join(dplyr::select(inner_join(dplyr::select(full_data_F5F6, ID, Final_weight), f5f6_snped, by = "ID"),
                           ID, pID ),
             t.df, by = "pID") %>%
  dplyr::select(ID, df$ID) %>% rename(id = ID)


write_csv(f5f6_solar, "./data/f5f6_genotypes_solar_chr6.csv")
