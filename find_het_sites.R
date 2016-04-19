source("./find_private_aleles.R")

by_line = laply(strsplit(colnames(snp_array), "_"), `[`, 1)
get_heterozygote = function(x) tapply(x, by_line, function(row) length(unique(row))) - 1

# is_heterozygote = tbl_df(ldply(snp_array_list,
#                                function(snp_array) aaply(snp_array, 1, get_heterozygote),
#                                .parallel = TRUE))[,-1]
# names(is_heterozygote) <- paste(unique(by_line), "hets", sep = "_")
# save(is_heterozygote, file = "./data/heterozygotes_position.Rdata")
load("./data/heterozygotes_position.Rdata")

just_snps <- bind_cols(just_snps, is_heterozygote)

just_hets <- just_snps %>% select(CHROM, POS, contains("_hets"))
just_hets

library(RcppRoll)

# het_array_list = dlply(just_hets, .(CHROM), function(x) select(x, contains("_hets")))
# mean_het = tbl_df(ldply(het_array_list,
#                  function(x) roll_mean(as.matrix(x),
#                                        floor(dim(x)[1]/10), fill = NA, align = "center"),
#                  .parallel = TRUE))
# mean_het = bind_cols(select(just_hets, CHROM, POS),
#                      mean_het[-1])
# names(mean_het)[3:8] <- unique(by_line)
# save(mean_het, file = "./data/mean_heterozygocity.Rdata")
load("./data/mean_heterozygocity.Rdata")

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