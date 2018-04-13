source("read_genotypes.R")

f5f6_snped = inner_join(full_data_F5F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID, contains("growth"), Final_weight)

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID, contains("growth"), Final_weight)

phenotypes = select(f6_snped, ID, Sex, Final_weight, growth_D3D7:growth_D21D28) %>% distinct(ID, .keep_all = TRUE)

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

idcf = cic(pedigree, phenotypes$ID, df = 1, ask = TRUE, msg = TRUE)
save(idcf, file = "./data/qtlrel/idcf_F6.Rdata")
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
               test = "None")

lrt$p
plot(lrt,gmap=map,main="Body Weight")
plotit(lrt$p)

plot(sort(lrt$p)~stats::ppoints(nrow(gwas)))
abline(0, 1)

idx<- match(colnames(gdatTmp),map$snp)
Tmp<- data.frame(chr=map$chr[idx],
                 dist=map$dist[idx],
                 y=lrt$p)
Tmp$chr <- reorder(Tmp$chr)
Tmp <- Tmp[order(Tmp$chr,Tmp$dist),] # order by chromosome and distance

plotit(Tmp, cv=12, main="Mapping Plot of Body Weight", xlab="Chromosome",
       ylab="LRT", col=as.integer(Tmp$ch)%%2+2,type="p",lty=2)

null_sims = nullSim(y=pdatTmp$Final_weight,
                    x=pdatTmp$Sex,
                    gdat=gdatTmp,
                    gmap = map,
                    vc = vc,
                    method = "permutation",
                    ntimes = 2)
null_sims
