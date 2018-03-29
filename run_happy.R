source("read_genotypes.R")

#detach("package:happy", unload = TRUE)
library(happy)

if(!require(kinship2)){install.packages("kinship2"); library(kinship2)}

.n = function(x) as.numeric(as.factor(x))
pedigree = as.data.frame(read.csv("./data/Intercross_pedigree2.csv")) %>% dplyr::rename(id = animal) %>% orderPed
ped2 = (left_join(dplyr::rename(pedigree, ID = id), dplyr::select(full_data, ID, Sex), by = "ID"))
for(i in 1:nrow(ped2)){
  if(is.na(ped2$Sex[i])){
    if(ped2$ID[i] %in% ped2$dam) ped2$Sex[i] = "F"
    else if(ped2$ID[i] %in% ped2$sire) ped2$Sex[i] = "M"
    else ped2$Sex[i] = "M"
  }
}
pedAll <- pedigree(id=ped2$ID,
                   dadid=ped2$sire, momid=ped2$dam,
                   sex=ped2$Sex)
A = Matrix(2*kinship(pedAll), sparse = TRUE)

h <- happy('./data/f5f6_genotypes_happy_chr6.PED', 'data/markers_strain_split_chr6.txt',
           generations=6,
           file.format = "ped",
           phase ="estimate" )

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
