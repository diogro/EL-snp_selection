library(readr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(doMC)
registerDoMC(4)

selected_snps = read_csv("./data/selected_snps.csv")

line_order = c("A13", "A31", "A41", "A23", "A22", "A42")[6:1]

sel_usnp_count <- selected_snps %>% count(snp_type, CHROM)
sel_usnp_count %>% spread(snp_type, n) %>% print(n = 21)
sel_usnp_quality <- selected_snps %>% group_by(CHROM) %>% summarize(min(QUAL), median(QUAL))

u_snps_plot_data <-
  selected_snps %>%
  select(CHROM, POS, snp_type, class) %>%
  mutate(POS = POS/1e6)

snp_order = as.character(unique(sel_usnp_count$snp_type))
snp_order = snp_order[order(nchar(snp_order))]
u_snps_plot_data$snp_type = factor(u_snps_plot_data$snp_type, levels = snp_order)
psnps_plots = llply(unique(selected_snps$CHROM),
                    function(current_chr)
                    {
                      counts = filter(sel_usnp_count, CHROM == current_chr)
                      u_line_order = as.character(unique(filter(sel_usnp_count, CHROM == current_chr)$snp_type))
                      u_line_order = u_line_order[order(nchar(u_line_order))]
                      counts = counts$n[match(u_line_order, counts$snp_type)]
                      legend = paste(u_line_order, "N =", counts)
                      private_alele_plot <-
                        u_snps_plot_data %>%
                        filter(CHROM == current_chr) %>%
                        ggplot(aes(snp_type, POS, color = snp_type)) +
                        geom_point(size = 0.2) + geom_point(size = 0.2, aes(0.5, POS)) +
                        coord_flip() + labs(x = "Line", y = "Chromossomal Position (Mb)") +
                        scale_color_discrete(labels = legend, name = "") +
                        ggtitle(current_chr)
                      save_plot(paste0("./data/jpegs/selected_usnps_", current_chr, ".png"),
                                private_alele_plot,
                                base_height = 6.5, base_aspect_ratio = 2.3)
                      return(private_alele_plot)
                    }, .parallel = TRUE)


current_chr = "chr1"
current_type = "A22"
counts = filter(sel_usnp_count, CHROM == current_chr)
u_line_order = as.character(unique(filter(sel_usnp_count, CHROM == current_chr)$snp_type))
u_line_order = u_line_order[order(nchar(u_line_order))]
counts = counts$n[match(u_line_order, counts$snp_type)]
legend = paste(u_line_order, "N =", counts)
private_alele_plot <-
  u_snps_plot_data %>%
  filter(CHROM == current_chr) %>%
  filter(grepl(current_type, snp_type)) %>%
  filter(class != "50_50") %>%
  mutate(n_class = as.numeric(factor(class)) + 2) %>%
  ggplot(aes(n_class, POS, color = snp_type)) +
  geom_point(size = 0.2) + geom_point(size = 0.2, aes(6, POS)) +
  coord_polar(theta = "y") + labs(x = "", y = "Chromossomal Position (Mb)") +
  scale_color_discrete(labels = legend, name = "") +
  scale_y_continuous(breaks = seq(0, 200, 10)) + scale_x_continuous(breaks = NULL, lim = c(0, 7)) +
  theme_bw() +
  ggtitle(current_type)
