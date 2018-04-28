source("read_genotypes.R")

#detach("package:happy", unload = TRUE)
library(happy)

f5f6_h <- happy('./data/happy_f5f6_genotypes.PED', 'data/happy_markers_strain.txt',
           generations=10,
           file.format = "ped",
           phase ="estimate" )
f6_h <- happy('./data/happy_f6_genotypes.PED', 'data/happy_markers_strain.txt',
           generations=10,
           file.format = "ped",
           phase ="estimate" )
gwas_happy_f5f6 = hfit(f5f6_h)
gwas_happy_f6 = hfit(f6_h)
save(f5f6_h, f6_h, gwas_happy_f5f6, gwas_happy_f6, file = "./data/HAPPY_f5f6_f6.Rdata")

gwas_happy_f5f6$model
gwas_table_happy = tbl_df(gwas_happy_f5f6$table)
gwas_table_happy$logP = as.numeric(gwas_table_happy$`additive logP`)
arrange(gwas_table_happy, logP)
names(gwas_table_happy)[2] = "rs"
gwh = inner_join(gwas_rsnp, select(gwas_table_happy, rs, logP), by = "rs")
gwh[is.na(gwh)] = 0
gwh$p = exp(-gwh$logP)



plot(gwas_happy)
happy::happyplot(gwas_happy)
marker = hdesign( h, h$markers[1000], model='additive')
x = tbl_df(marker) %>%
  mutate(ID = as.numeric(h$subjects))
names(full_data)
data = left_join(x, full_data_F5F6, by = "ID") %>%
  mutate(animal = ID) %>%
  filter(Gen == "F6") %>%
  dplyr::select(ID, animal, Final_weight, Sex, Litter_size_birth, Birth_litter_size_weaning, Foster_litter_size_weaning, line_order) %>%
  filter(!is.na(Final_weight))

fixed_term = "Sex + Litter_size_birth + Birth_litter_size_weaning + Foster_litter_size_weaning"
random_term = "(1|ID)"
marker_term = paste(line_order, collapse = " + ")
null_formula = paste0("Final_weight ~ ", fixed_term)
null_random_formula = paste(null_formula, random_term, sep = " + ")
marker_formula = paste(null_formula, marker_term, sep = " + ")
marker_random_formula = paste(null_random_formula, marker_term, sep = " + ")

ped2$animal <- as.factor(ped2$ID)
prior_uni <- list(G = list(G1 = list(V = 1, nu = 0.002)),
                  R = list(V = 1, nu = 0.002))
null_mcmc = MCMCglmm(as.formula(null_formula),
                     random = ~ animal,
                     prior = prior_uni,
                     pedigree = ped2[,c("animal", "dam", "sire")],
                     data = as.data.frame(data))
marker_mcmc = MCMCglmm(as.formula(marker_formula),
                       random = ~ animal,
                       prior = prior_uni,
                       pedigree = ped2[,c("animal", "dam", "sire")],
                       data = as.data.frame(data))
summary(null_mcmc)
summary(marker_mcmc)

var(lm(as.formula(null_formula), data = data) %>% residuals)
null_lm = relmatLmer(as.formula(null_random_formula),   data = data, relmat = list(ID = A), REML=F,
                     control = lmerControl(calc.derivs = FALSE))
mark_lm = relmatLmer(as.formula(marker_random_formula), data = data, relmat = list(ID = A), REML=F,
                     control = lmerControl(calc.derivs = FALSE))
update_formula = paste0(". ~ . + ", marker_term)
mark_lm = update(null_lm, as.formula(update_formula), control = lmerControl(calc.derivs = FALSE))


summary(null_lm)
summary(mark_lm)

relgrad <- with(mark_lm@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))

lmerTest::anova(null_lm, mark_lm)
?calcSatterth
