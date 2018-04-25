source("./read_phenotypes.R")

raw_gen_missing = read_table2("./data/plink_files/atchley_missing_text.vcf", comment = "##") %>% arrange(`#CHROM`, POS)
raw_gen_imputed = read_table2("./data/plink_files/atchley_imputed_text.vcf", comment = "##") %>% arrange(`#CHROM`, POS)
raw_gen_missing = raw_gen_missing[match(raw_gen_imputed$ID, raw_gen_missing$ID),]

col_ID = colnames(raw_gen_imputed)[-c(1:9)]

col = 1000
diff_list = llply(seq_along(col_ID), function(col){
  diff = which(!raw_gen_imputed[,col_ID[col]] == raw_gen_missing[,col_ID[col]])
  data.frame(raw_gen_imputed[diff,1:3],
             imp = raw_gen_imputed[diff,col_ID[col]],
             mis = raw_gen_missing[diff,col_ID[col]])
})

ndiffs = laply(diff_list, nrow)
summary(ndiffs)
hist(ndiffs, breaks = 200)
which(ndiffs == median(ndiffs))
which(ndiffs == 201)
diff_list[813]
ggplot(raw_gen_imputed, aes(POS, `#CHROM`)) + geom_point(size = 0.01) +
  geom_point(data = diff_list[[1059]], aes(POS, X.CHROM), color = "tomato3", size = 1)


sum(ndiffs > 200)/length(ndiffs)
