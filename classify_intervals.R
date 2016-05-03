source("./find_het_sites.R")

classifyInterval <- function(snp_density) {
  snp_density = bind_cols(select(snp_density_binary, CHROM, POS), snp_density[-1])
  names(snp_density) = names(snp_density_binary)
  snp_density = select_(snp_density, "CHROM","POS",
                       "A22_hets", "A22_het_p_A23", "A22_het_p_A42", "A22_pri_het",
                       "A22", "A22A42", "A22A23")
  snp_density = na.omit(snp_density)

  mask_het = c(FALSE, FALSE, FALSE, rep(TRUE, 6))
  mask_non_het = c(FALSE, FALSE, FALSE, rep(FALSE, 3), rep(TRUE, 3))
  classify_A22 = function(x){
    id = x[,1:2]
    if(x$A22_hets < 0.3){
      dx = sort(x[,mask_non_het], decreasing = TRUE)
      out = bind_cols(id, data_frame(is_het = 0, alele_1 = names(dx)[1], alele_2 = NA))
    }  else{
      dx = sort(x[,mask_het], decreasing = TRUE)
      out = bind_cols(id, data_frame(is_het = 1, alele_1 = names(dx)[1], alele_2 = names(dx)[2]))
    }
    if(all(dx < 0.001)){out$alele_1 = NA; out$alele_2 = NA}
    return(out)
  }
  snp_density_list = dlply(snp_density, .(CHROM), identity)
  interval_class = ldply(snp_density_list,
                         function(x) ddply(x, .(POS), classify_A22, .progress = "text"),
                         .parallel = TRUE)
  interval_class
}
interval_class_100 = classifyInterval(snp_density_100)
interval_class_1000 = classifyInterval(snp_density_1000)
#save(interval_class_100, interval_class_1000, file = "data/intervalClass.Rdata")
#load("./data/intervalClass.Rdata")

interval_class_plots_100 = llply(unique(just_snps$CHROM),
                             function(current_chr)
                             {
                               interval_class_plot <-
                                 interval_class_100 %>%
                                 filter(CHROM == current_chr) %>%
                                 ggplot() +
                                 geom_point(size = 0.3, alpha = 0.1, aes(alele_1,
                                                                         POS,
                                                                         color = alele_1)) +
                                 geom_point(size = 0.3, alpha = 0.1, aes(alele_2,
                                                                         POS,
                                                                         color = alele_2)) +
                                 coord_flip() + labs(x = "Class",
                                                     y = "Chromossomal Position (Mb)") +
                                 ggtitle(current_chr)
                               save_plot(paste0("./data/jpegs/interval_class_",
                                                current_chr, "_res100.png"),
                                         interval_class_plot,
                                         base_height = 5, base_aspect_ratio = 2)
                               return(interval_class_plot)
                             }, .parallel = TRUE)

interval_class_plots_1000 = llply(unique(just_snps$CHROM),
                             function(current_chr)
                             {
                               interval_class_plot <-
                                 interval_class_1000 %>%
                                 filter(CHROM == current_chr) %>%
                                 ggplot() +
                                 geom_point(size = 0.3, alpha = 0.1, aes(alele_1,
                                                                         POS,
                                                                         color = alele_1)) +
                                 geom_point(size = 0.3, alpha = 0.1, aes(alele_2,
                                                                         POS,
                                                                         color = alele_2)) +
                                 coord_flip() + labs(x = "Class",
                                                     y = "Chromossomal Position (Mb)") +
                                 ggtitle(current_chr)
                               save_plot(paste0("./data/jpegs/interval_class_",
                                                current_chr, "_res1000.png"),
                                         interval_class_plot,
                                         base_height = 5, base_aspect_ratio = 2)
                               return(interval_class_plot)
                             }, .parallel = TRUE)
