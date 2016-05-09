source("./find_private_aleles.R")

by_line = laply(strsplit(colnames(snp_array_list[[1]]), "_"), `[`, 1)
get_heterozygote = function(x) tapply(x, by_line, function(row) length(unique(row))) - 1

is_heterozygote = tbl_df(ldply(snp_array_list,
                               function(snp_array) aaply(snp_array, 1, get_heterozygote),
                               .parallel = TRUE))[,-1]
names(is_heterozygote) <- paste(unique(by_line), "hets", sep = "_")
save(is_heterozygote, file = "./data/heterozygotes_position_3_strains.Rdata")
load("./data/heterozygotes_position_3_strains.Rdata")

just_snps <- bind_cols(just_snps, is_heterozygote)

just_hets <- just_snps %>% select(CHROM, POS, contains("_hets"))
just_hets

library(RcppRoll)

het_array_list = dlply(just_hets, .(CHROM), function(x) select(x, contains("_hets")))
mean_het = tbl_df(ldply(het_array_list,
                 function(x) roll_mean(as.matrix(x),
                                       floor(dim(x)[1]/1000), fill = NA, align = "center"),
                 .parallel = TRUE))
mean_het = bind_cols(select(just_hets, CHROM, POS),
                     mean_het[-1])
names(mean_het)[3:6] <- unique(by_line)
save(mean_het, file = "./data/mean_heterozygocity_3_strains.Rdata")
load("./data/mean_heterozygocity_3_strains.Rdata")

m_mean_het = melt(mean_het, id.vars = c("CHROM", "POS"))
m_mean_het = mutate(m_mean_het,
                    variable = factor(variable, levels = line_order),
                    POS = POS/1e6)

het_plots = llply(unique(just_snps$CHROM),
                  function(current_chr)
                  {
                    het_plot = ggplot(filter(m_mean_het, CHROM == current_chr),
                                      aes(POS, value, group = variable, color = variable)) +
                      geom_line() + labs(y = "Mean heterozygocity",
                                         x = "Chromossomal Position (Mb)") +
                      scale_color_discrete(name = "") + ylim(0, 1) +
                      ggtitle(current_chr)
                    save_plot(paste0("./data/jpegs/het_", current_chr, ".png"),
                              het_plot,
                              base_height = 5, base_aspect_ratio = 2)
                    return(het_plot)
                  }, .parallel = TRUE)

snp_density_binary <- select(just_hets, CHROM, POS, A22_hets)
snp_types = as.character(unique(just_snps$pu_line)[-1])
for(current in snp_types){
    new_var = paste0(current)
    new_col = as.numeric(just_snps$pu_line == current)
    new_col[is.na(new_col)] = 0
    snp_density_binary[[new_var]] = new_col
}

calculate_snpDensity = function (resolution){
  snp_density_list = dlply(snp_density_binary, .(CHROM), function(x) select(x, -CHROM, -POS))
  snp_density = tbl_df(ldply(snp_density_list,
                             function(x) roll_mean(as.matrix(x),
                                                   floor(dim(x)[1])/resolution, fill = NA, align = "center"),
                             .parallel = TRUE))
  snp_density = bind_cols(select(snp_density_binary, CHROM, POS), snp_density[-1])
  names(snp_density) = names(snp_density_binary)
  snp_density
}
snp_density_10   = calculate_snpDensity(10)
snp_density_100  = calculate_snpDensity(100)
snp_density_1000 = calculate_snpDensity(1000)
save(snp_density_10, snp_density_100, snp_density_1000, file = "data/snp_density_3_strains.Rdata")
load("data/snp_density_3_strains.Rdata")

snp_density =  (snp_density_10)
create_snpDesityPlot <- function(resolution, snp_density) {
  snp_density = select(snp_density, -A13)

  m_snp_density = melt(snp_density,
                       id.vars = c("CHROM", "POS"))
  m_snp_density = mutate(m_snp_density,
                      POS = POS/1e6)
  m_snp_density[["class"]] <- c(NA, "A22r")[grepl("A22", m_snp_density$variable)+1]
  m_snp_density[["class"]][m_snp_density$variable == "A42"] <- "A42"
  m_snp_density[["class"]][m_snp_density$variable == "A23"] <- "A23"
  m_snp_density[["class"]][m_snp_density$variable == "A22_hets"] <- "A22_hets"

  levels(m_snp_density$variable)
  class_colors = c("black", "red", "green", "blue", "darkblue", "darkred", "black", "red", "blue")
  names(class_colors) <- levels(m_snp_density$variable)
  current_chr = "chr3"
  snp_density_plots = llply(unique(just_snps$CHROM),
                    function(current_chr)
                    {
                        snp_density_plot <-
                          ggplot(filter(m_snp_density, CHROM == current_chr),
                                          aes(POS, log1p(value), group = variable, color = variable)) +
                    geom_line() + labs(y = "Mean density",
                                       x = "Chromossomal Position (Mb)") +
                    scale_color_manual(name = "", values = class_colors) + #ylim(0, 0.001) +
                    facet_wrap(~class, nrow = 4, scales = "free_y") +
                    ggtitle(current_chr)
                save_plot(paste0("./data/jpegs/snp_density", current_chr, "res_", resolution, ".png"),
                          snp_density_plot,
                          base_height = 8, base_aspect_ratio = 2)
                return(snp_density_plot)
                    }, .parallel = TRUE)
}
create_snpDesityPlot(10  , snp_density_10)
create_snpDesityPlot(100 , snp_density_100)
create_snpDesityPlot(1000, snp_density_1000)
