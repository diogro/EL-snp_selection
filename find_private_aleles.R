source("./read_raw_snp_data.R")

is_private = function(x){
  tx = table(x)
  if(any(tx == 2)){
    rare = names(tx)[which(min(tx) == tx)]
    lines = names(x)[x == rare]
    line = unique(laply(strsplit(lines, "_"), `[`, 1))
    if(length(line) == 1) return(c(TRUE, line))
    else return(c(FALSE, NA))
  }
  else return(c(FALSE, NA))
}

by_line = laply(strsplit(colnames(snp_array), "_"), `[`, 1)
is_private2 = function(x){
  tapply(x, by_line, unique)
  #TODO
}

#get_private = aaply(snp_array, 1, is_private, .parallel = TRUE)
#save(get_private, file = "./data/private_snp_position.Rdata")
load("./data/private_snp_position.Rdata")

just_snps = mutate(just_snps,
                   is_private = as.logical(get_private[,1]),
                   p_line = factor(get_private[,2], levels = line_order))

snp_count <-
  just_snps %>%
  filter(is_private) %>%
  count(p_line, CHROM)

# snps_plot <-
#   just_snps %>%
#   filter(is_private) %>%
#   select(CHROM, POS, p_line) %>%
#   mutate(POS = POS/1e6)
# for(current_chr in unique(just_snps$CHROM)){
#   counts = filter(snp_count, CHROM == current_chr)
#   counts = counts$n[which(line_order == filter(snp_count, CHROM == current_chr)$p_line)]
#   legend = paste(line_order, "N =", counts)
#   private_alele_plot <-
#     snps_plot %>%
#     filter(CHROM == current_chr) %>%
#     ggplot(aes(p_line, POS, color = p_line)) +
#     geom_point(size = 0.3) + geom_point(size = 0.2, aes(0.5, POS)) +
#     coord_flip() + labs(x = "Line", y = "Chromossomal Position (Mb)") +
#     scale_color_discrete(labels = legend, name = "") +
#     ggtitle(current_chr)
#    save_plot(paste("./data/jpegs/new_", current_chr, ".png"),
#              private_alele_plot,
#              base_height = 5, base_aspect_ratio = 2)
# }