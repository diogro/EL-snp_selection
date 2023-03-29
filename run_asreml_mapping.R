source("read_genotypes.R")
library(asreml)
library(corrplot)

f5f6_snped = inner_join(full_data_F5F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID, contains("growth"), Final_weight)

growthF5F6 = select(f5f6_snped, ID, Sex, Final_weight, growth_D0D14:growth_D42D56) %>% distinct(ID, .keep_all = TRUE)

growthF5F6 = growthF5F6[complete.cases(growthF5F6[,growth_traits]),]

pedigree = as.data.frame(read.csv("./data/Intercross_pedigree2.csv")) %>%
  rename(id = animal) %>% orderPed
pedigree$id = factor(pedigree$id)
pedigree$dam = factor(pedigree$dam)
pedigree$sire = factor(pedigree$sire)
names(pedigree)[1] = "ID"

ginverse =  asreml.Ainverse(pedigree)
Ainv = asreml.Ainverse(pedigree)$ginv

n_traits = length(growth_traits)

growthF5F6_std = growthF5F6
growthF5F6_std[,growth_traits] = scale(growthF5F6_std[,growth_traits])
growthF5F6_sd = apply(growthF5F6[,growth_traits], 2, sd)

data_growth = as.data.frame(growthF5F6_std)
names(data_growth)

g.sv = asreml(cbind(growth_D0D14, growth_D14D28, growth_D28D42, growth_D42D56) ~ trait + trait:Sex,
                           random = ~ us(trait):ped(ID),
                           rcov = ~ units:us(trait),
                           ginverse = list(ID = Ainv),
                           data = data_growth, start.values = TRUE)
sv = g.sv$gammas.table
GR = read.csv("mixed_model_start.csv")
sv$Value[1:10] = GR$G
sv$Value[12:21] = GR$R
growth_animal_asr = asreml(cbind(growth_D0D14, growth_D14D28, growth_D28D42, growth_D42D56) ~ trait + trait:Sex,
                           random = ~ us(trait):ped(ID),
                           rcov = ~ units:diag(trait),
                           ginverse = list(ID = Ainv),
                           data = data_growth,
                           G.param = sv)

G = matrix(NA, n_traits, n_traits)
G[upper.tri(G, diag = TRUE)] = growth_animal_asr$G.param$`trait:ped(ID)`$trait$initial
G[lower.tri(G)] = t(G)[lower.tri(G)]
G = G * outer(growthF5F6_sd, growthF5F6_sd)
corrG = cov2cor(G)
colnames(corrG) = c("0 to 14\ndays", "14 to 28\ndays", "28 to 42\ndays", "42 to 56\ndays")
diag(corrG) = 0
corrplot.mixed(corrG, upper = "ellipse", mar = c(0, 0, 0, 0), cl.lim = c(-0.8, 0.8), addgrid.col = "black", is.corr = FALSE)





if(!require(data.table)){install.packages("data.table"); library(data.table)}

ped3 = ped2
ped3$pedID = ped3$ID
is.na(ped3)
ped3[is.na(ped3)] = 0
ped3 = select(ped3, pedID, ID, sire, dam)
write_delim(ped3, "./data/Intercross_pedigree3.csv", delim = " ", col_names = FALSE)

dt_gen = data.table(gen[,-c(1:4)])
for(col in names(dt_gen)) set(dt_gen, i=which(dt_gen[[col]]=="1/1"), j=col, value=1)
for(col in names(dt_gen)) set(dt_gen, i=which(dt_gen[[col]]=="0/1"), j=col, value=2)
for(col in names(dt_gen)) set(dt_gen, i=which(dt_gen[[col]]=="0/0"), j=col, value=3)
dt_gen <- dt_gen[, lapply(.SD, as.numeric)]
qtl_rel_gen = t(dt_gen)
dimnames(qtl_rel_gen)[[2]] = gen$ID

dim(qtl_rel_gen)
gdatF6 = as.matrix(qtl_rel_gen)

map = gen %>%
  select(ID, chr, gpos, pos) %>%
  rename(snp = ID,
         dist = gpos,
         phyPos = pos)

#idcf = cic(pedigree, phenotypes$ID, df = 1, ask = TRUE, msg = TRUE)
#save(idcf, file = "./data/qtlrel/idcf_F6.Rdata")
load("./data/qtlrel/idcf_F6.Rdata")
gmF6 <- genMatrix(idcf)

idx <- !is.na(phenotypes$Final_weight)
pdatTmp <- phenotypes[idx,]
gdatTmp <- gdatF6[match(pdatTmp$ID,rownames(gdatF6)),]
ii<- match(pdatTmp$ID,rownames(gmF6$AA))
vc = estVC(pdatTmp$Final_weight,pdatTmp$Sex,
      v=list(AA=gmF6$AA[ii,ii],
             DD=gmF6$DD[ii,ii],
             HH=NULL,
             AD=NULL,
             MH=NULL,
             E=diag(nrow(pdatTmp))))

lrt <- scanOne(y=pdatTmp$Final_weight,
               x=pdatTmp$Sex,
               gdat=gdatTmp,
               vc=vc,
               test = "Chisq")

plot(-log10(lrt$p))
plot(-log10(gwas_rped$p_lrt))
plot(lrt,gmap=map,main="Body Weight")
plotit(lrt$p)

plot((lrt$p)~(gwas_rped$p_lrt))
plot(obs_rped, obs_qtl_rel)
abline(0, 1)

plot(sort(lrt$p)~stats::ppoints(nrow(gwas)))
abline(0, 1)

idx<- match(colnames(gdatTmp),map$snp)
Tmp<- data.frame(chr=map$chr[idx],
                 dist=map$dist[idx],
                 y=lrt$p)
Tmp$chr <- reorder(Tmp$chr)
Tmp <- Tmp[order(Tmp$chr,Tmp$dist),] # order by chromosome and distance

plotit(Tmp, cv=5, main="Mapping Plot of Body Weight", xlab="Chromosome",
       ylab="LRT", col=as.integer(Tmp$ch)%%2+2,type="p",lty=2)

null_sims = nullSim(y=pdatTmp$Final_weight,
                    x=pdatTmp$Sex,
                    gdat=gdatTmp,
                    gmap = map,
                    vc = vc,
                    method = "permutation",
                    ntimes = 2)
null_sims
