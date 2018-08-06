source("./read_phenotypes.R")
source("./modular_classes.R")

library(evolqg)
#install.packages("superheat")
library(superheat)
vectorCor = function(x, y) Normalize(x) %*% Normalize(y)

gwas_rloco = ldply(1:20, function(i) read_tsv(paste0("./output/growth_r-loco_multivariate_2w_t1-4_chr",i,".assoc.txt")))
gwas_rsnp = ldply(1:20, function(i) read_tsv(paste0("./output/growth_r-snp_multivariate_2w_t1-4_chr",i,".assoc.txt")))
gwas_rped = ldply(1:20, function(i) read_tsv(paste0("./output/growth_r-ped_multivariate_2w_t1-4_chr",i,".assoc.txt")))

plot.inflation <- function (x, size = 2) {

  # Get the number of p-values.
  n <- length(x)

  # Compute the negative log10(p-values), and sort them from largest
  # to smallest.
  y <- rev(sort(-log10(x)))

  # Create the q-q plot.
  return(ggplot(data.frame(x = -log10((1:n)/n),y = y),aes(x = x,y = y)) +
           geom_abline(intercept = 0,slope = 1,color = "magenta") +
           geom_point(color = "dodgerblue",shape = 20,size = 2) +
           labs(x = "Expected -log10 p-value",y = "Observed -log10 p-value") +
           theme(axis.line = element_blank()))
}

plot.inflation(gwas_rped$p_wald)
plot.inflation(gwas_rloco$p_wald)
plot.inflation(gwas_rsnp$p_wald)

table(gwas_rsnp$p_wald < 5.17E-7)
gwas_rped[which(gwas_rped$p_wald < 5.17E-7),]
(gwas_growth_p_ped = ggman(gwas_rped, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_wald", relative.positions = TRUE, title = "GEMMA ped growth", sigLine = -log10(2.6e-5), pointSize = 1))
(gwas_growth_p_snp = ggman(gwas_rsnp, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_wald", relative.positions = TRUE, title = "GEMMA snp growth", sigLine = -log10(2.6e-5), pointSize = 1))
(gwas_growth_p_loco = ggman(gwas_rloco, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_wald", relative.positions = TRUE, title = "GEMMA snp growth", sigLine = -log10(2.6e-5), pointSize = 1))


positions = read.table("./data/marker_positions.txt", header = FALSE, stringsAsFactors = FALSE) %>% tbl_df
names(positions) = c("chr2", "ps", "chr", "gpos")
positions$chr2 = NULL
positions$chr = as.integer(positions$chr)

gwas = gwas_rped
gwas = inner_join(gwas, positions, by = c("chr","ps")) %>%
  filter(p_wald < genome_wide_p)

getChoosenSNPs <- function(gwas) {
  gwas = inner_join(gwas, positions, by = c("chr","ps")) %>%
    filter(p_wald < genome_wide_p)
  choosen_snp_list = list()
  k = 1
  for(c in unique(gwas$chr)){
    while(TRUE){
      min_p = which.min(gwas[gwas$chr == c, "p_wald"])
      choosen_snp = as.character(gwas[gwas$chr == c,"rs"][min_p])
      if(length(choosen_snp) > 0){
        choosen_snp_list[[k]] = choosen_snp
        k = k + 1
      } else {
        break
      }
      choosen_pos = gwas[gwas$rs == choosen_snp, "gpos"]
      removable_snps = gwas[gwas$chr == c, "rs"][abs(gwas[gwas$chr == c, "gpos"] - choosen_pos) < 20]
      gwas = filter(gwas, !rs %in% removable_snps)
    }
  }
  unlist(choosen_snp_list)
}
genome_wide_p = 1e-4
g_rloco_sig = gwas_rloco %>% filter(rs %in% getChoosenSNPs(.))

p1 = ggman(gwas_rloco, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_wald",
      relative.positions = TRUE, title = "GEMMA rLOCO", sigLine = -log10(genome_wide_p),
      pointSize = 1) 
ggmanLabel(p1, labelDfm = g_rloco_sig, snp = "rs", label = "rs")

ad = as.matrix(g_rloco_sig[,paste("beta", 1:4, sep = "_")])
ad_classes = apply(ad, 1, classifyVector)
rownames(ad) = ad_classes
colnames(ad) = paste(c("0 to 14", "14 to 28", "28 to 42", "42 to 56"), "\n days")

png("~/Dropbox/labbio/posters/2018\ -\ 07\ -\ 19\ -\ Evolution2018/additive_effects_bars.png", width = 700, height = 900)
superheat(10*ad, 
          order.rows = order(ad_classes, decreasing = TRUE),
          heat.pal = c("#b35806", "white", "#542788"),
          heat.pal.values = c(0, 0.47, 1),
          X.text = round(as.matrix(10*ad), 1),
          X.text.size = 4, 
          legend = FALSE,
          grid.vline = FALSE,
          membership.rows = rownames(ad),
          yr = apply(10*ad, 1, Norm),
          yr.plot.type = "bar",
          yr.axis.name = "Vector size",
          yr.bar.col = "black",
          yr.obs.col = rep("grey", nrow(ad)))
dev.off()
png("~/Dropbox/labbio/posters/2018\ -\ 07\ -\ 19\ -\ Evolution2018/additive_effects_points.png", width = 1550, height = 1800)
superheat(10*ad, 
          order.rows = order(ad_classes, decreasing = TRUE),
          heat.pal = c("#b35806", "white", "#542788"),
          heat.pal.values = c(0, 0.47, 1),
          X.text = round(as.matrix(10*ad), 1),
          X.text.size = 15, 
          legend = FALSE,
          grid.vline = FALSE,
          membership.rows = rownames(ad),
          left.label.text.size = 12, 
          bottom.label.text.size = 15,
          yr = apply(10*ad, 1, Norm),
          yr.line.col = "black",
          yr.line.size = 10,
          yr.axis.name = "Vector\nsize",
          yr.obs.col = rep("black", nrow(ad)),
          yr.point.size = 8,
          yr.axis.size = 20,
          yr.axis.name.size = 40)
dev.off()
