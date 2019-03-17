source("./read_phenotypes.R")
source("./modular_classes.R")

library(evolqg)
if(!require(superheat)){install.packages("superheat"); library(superheat)}
if(!require(corrplot)){install.packages("corrplot"); library(corrplot)}
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

per_chrom_tresholds_files = paste0("./data/plink_files/per_chrom/gec_chr", 1:20, ".sum")
per_chrom_p = data.frame(chr = 1:20, treshold = ldply(per_chrom_tresholds_files, read_tsv)$Significant_P_Value)

gwas = gwas_rloco
gwas = inner_join(gwas, positions, by = c("chr","ps"))
gwas = ldply(1:20, function(i) filter(gwas, chr == i) %>% filter(p_wald < per_chrom_p[i,2]))

getChoosenSNPs <- function(gwas) {
  gwas = inner_join(gwas, positions, by = c("chr","ps"))
  gwas = ldply(1:20, function(i) filter(gwas, chr == i) %>% filter(p_wald < per_chrom_p[i,2]))
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
genome_wide_p = 0.1/5388.44

g_rloco_sig = gwas_rloco %>% filter(rs %in% getChoosenSNPs(.))

p1 = ggman(gwas_rloco, snp = "rs", bp = "ps", chrom = "chr", pvalue = "p_wald",
      relative.positions = TRUE, title = "", sigLine = -log10(genome_wide_p),
      pointSize = 1)
manhattan = ggmanLabel(p1, labelDfm = g_rloco_sig, snp = "rs", label = "rs") + theme_cowplot()
ggmanHighlight(p1, highlight = toy.highlights)

save_plot("~/Dropbox/labbio/articles/TeseDoutorado/chapter_atchley/media/selected_snps.png", manhattan, base_height = 7, base_aspect_ratio = 2)

ad = as.matrix(g_rloco_sig[,paste("beta", 1:4, sep = "_")])
ad_classes = apply(ad, 1, classifyVector)
rownames(ad) = ad_classes
colnames(ad) = paste(c("0 to 14", "14 to 28", "28 to 42", "42 to 56"), "\n days")

library(xtable)
ad_table = cbind(left_join(g_rloco_sig[,c("chr", "rs", "ps", "p_wald")], positions, by = c("chr","ps")), ad_classes)
ad_table = ad_table[c("rs", "chr", "ps", "gpos", "ad_classes", "p_wald")]
ad_table$p_wald = -log10(ad_table$p_wald)
xtable(ad_table)

png("~/Dropbox/labbio/articles/TeseDoutorado/chapter_atchley/media/additive_effects_points.png", width = 1550, height = 1800)
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

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Foster_litter_size_weaning, contains("growth"), Final_weight)

phenotypes = select(f6_snped, ID, Sex, Final_weight,
                    Litter_size_birth,
                    Birth_litter_size_weaning,
                    Foster_litter_size_weaning,
                    growth_traits) %>% distinct(ID, .keep_all = TRUE)
trait_sd = apply(phenotypes[,growth_traits], 2, sd, na.rm = TRUE)

current_snp = ad_table$rs[1]
f6_ids = inner_join(full_data_F6, IDs, by = "ID")$ID
gen_df = as.data.frame(gen)
x = unlist(as.list(gen_df[gen_df$ID == current_snp, f6_ids]))
table(x)
names(x) = NULL
unique(x)
class(x[4])

a_effects_list = dlply(g_rloco_sig, .(rs), function(x) list(ID = x$rs,
                                           ad_effect = as.vector(x[,paste("beta", 1:4, sep = "_")])))

calcVa = function(current_snp){
  a_effects = as.numeric(a_effects_list[[current_snp]]$ad_effect)
  focal_marker_col = unlist(as.list(gen_df[gen_df$ID == current_snp, f6_ids]))
  n = length(focal_marker_col)
  genotype_freq = table(focal_marker_col)/n
  # Variance due to focal marker
  q = genotype_freq[1] + 1/2 * genotype_freq[2]
  p = genotype_freq[3] + 1/2 * genotype_freq[2]
  # additive contribution to Va
  V_a = 2*p*q * outer(a_effects, a_effects)
  V_a
}
Va_list = llply(ad_table$rs, calcVa)
llply(Va_list, cov2cor)

classes = levels(ad_table$ad_classes)
current_class = classes[[4]]
Va = list()
for(i in seq_along(classes))
  Va[[i]] = Reduce("+", Va_list[classes == classes[[i]]])
names(Va) = classes
Va[[6]] = Reduce("+", Va_list)

png("~/Dropbox/labbio/articles/TeseDoutorado/chapter_atchley/media/additive_matrices_by_class.png", width = 900, height = 1400)
layout(matrix(c(1, 2, 3, 4, 5, 5), 3, 2, byrow = TRUE))
corrplot.mixed(cov2cor(Va[["Modular"]]), upper = "ellipse", main = "Modular",number.cex=3,tl.cex=3, cex.main=3, mar = c(2, 2, 2, 2))
corrplot.mixed(cov2cor(Va[["Intra\nmodule"]]), upper = "ellipse", main = "Intra-module",number.cex=3,tl.cex=3,cex.main=3, mar = c(2, 2, 2, 2))
corrplot.mixed(cov2cor(Va[["Integrated"]]), upper = "ellipse", main = "Integrated", number.cex=3,tl.cex=3,cex.main=3, mar = c(2, 2, 2, 2))
corrplot.mixed(cov2cor(Va[["Antagonistic"]]), upper = "ellipse", main = "Antagonistic", number.cex=3,tl.cex=3,cex.main=3, mar = c(2, 2, 2, 2))
corrplot.mixed(cov2cor(Va[[6]]), upper = "ellipse", main = "Total additive genetic", number.cex=3,tl.cex=3,cex.main=3, mar = c(2, 2, 2, 2))
dev.off()
