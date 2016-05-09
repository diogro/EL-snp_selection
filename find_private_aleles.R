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

x = snp_array_list[[1]][1,]
A22_alele1_mask = as.character((names(x) == "A22_1") + 1)
A22_alele2_mask = (names(x) == "A22_2") + 1
get_A22_p_het = function(x){
  if(all(tapply(x, A22_alele1_mask, function(row) length(unique(row))) == c(1, 1)))
    return(TRUE)
  else if(all(tapply(x, A22_alele2_mask, function(row) length(unique(row))) == c(1, 1)))
    return(TRUE)
  else
    return(FALSE)
}

A22A23_mask = as.character((grepl("A22", names(x)) | grepl("A23", names(x))) + 1)
get_A22A23_p = function(x){
  if(all(tapply(x, A22A23_mask, function(row) length(unique(row))) == c(1, 1)))
    return(TRUE)
  else
    return(FALSE)
}

A22A42_mask = as.character((grepl("A22", names(x)) | grepl("A42", names(x))) + 1)
get_A22A42_p = function(x){
  if(all(tapply(x, A22A42_mask, function(row) length(unique(row))) == c(1, 1)))
    return(TRUE)
  else
    return(FALSE)
}

A22_alele1_A23_mask = as.character((grepl("A22_1", names(x)) | grepl("A23", names(x))) + 1)
A22_alele2_A23_mask = as.character((grepl("A22_2", names(x)) | grepl("A23", names(x))) + 1)
get_A22_het_p_23 = function(x){
  if(all(tapply(x, A22_alele1_A23_mask, function(row) length(unique(row))) == c(1, 1)))
    return(TRUE)
  else if(all(tapply(x, A22_alele2_A23_mask, function(row) length(unique(row))) == c(1, 1)))
    return(TRUE)
  else
    return(FALSE)
}

A22_alele1_A42_mask = as.character((grepl("A22_1", names(x)) | grepl("A42", names(x))) + 1)
A22_alele2_A42_mask = as.character((grepl("A22_2", names(x)) | grepl("A42", names(x))) + 1)
get_A22_het_p_42 = function(x){
  if(all(tapply(x, A22_alele1_A42_mask, function(row) length(unique(row))) == c(1, 1)))
    return(TRUE)
  else if(all(tapply(x, A22_alele2_A42_mask, function(row) length(unique(row))) == c(1, 1)))
    return(TRUE)
  else
    return(FALSE)
}

by_line = laply(strsplit(colnames(snp_array_list[[1]]), "_"), `[`, 1)
is_any_kind = function(x){
  tx = table(x)
  if(any(tx == 2)){
    rare = names(tx)[which(min(tx) == tx)]
    lines = names(x)[x == rare]
    line = unique(laply(strsplit(lines, "_"), `[`, 1))
    if(length(line) == 1) return(c(TRUE, line))
  } else if(any(tx == 1)){
    if(get_A22_p_het(x)) return(c(TRUE, "A22_pri_het"))
  } else if(any(tx == 4)){
    if(get_A22A42_p(x)) return(c(TRUE, "A22A42"))
    else if(get_A22A23_p(x)) return(c(TRUE, "A22A23"))
  } else if(any(tx == 3)){
    if(get_A22_het_p_42(x)) return(c(TRUE, "A22_het_p_A42"))
    else if(get_A22_het_p_23(x)) return(c(TRUE, "A22_het_p_A23"))
  }
  return(c(FALSE, NA))
}

get_private = tbl_df(ldply(snp_array_list,
                    function(snp_array) adply(snp_array, 1, is_private),
                    .parallel = TRUE))
save(get_private, file = "./data/private_snp_position.Rdata")
#load("./data/private_snp_position.Rdata")

is_usefull = tbl_df(ldply(snp_array_list,
                          function(snp_array) adply(snp_array, 1, is_any_kind),
                          .parallel = TRUE))
save(is_usefull, file = "./data/usefull_snp_position.Rdata")
#load("./data/usefull_snp_position.Rdata")

line_order = c("A13", "A31", "A41", "A23", "A22", "A42")[6:1]
just_snps = mutate(just_snps,
                   is_private = as.logical(get_private$V1),
                   p_line = factor(get_private$V2, levels = line_order),
                   is_usefull = as.logical(is_usefull$V1),
                   pu_line = as.factor(is_usefull$V2))

private_snps = just_snps %>% filter(is_usefull) %>% select(-is_usefull)

p_snp_count <-
  private_snps %>%
  filter(is_private) %>%
  count(p_line, CHROM) 
p_snp_count %>% spread(p_line, n) %>% print(n = 21)

#snps_plot_data <-
#  private_snps %>%
#  filter(is_private) %>%
#  select(CHROM, POS, p_line) %>%
#  mutate(POS = POS/1e6)
#psnps_plots = llply(unique(just_snps$CHROM),
#                    function(current_chr)
#                    {
#                      counts = filter(p_snp_count, CHROM == current_chr)
#                      counts = counts$n[which(line_order == filter(p_snp_count,
#                                                                   CHROM == current_chr)$p_line)]
#                      legend = paste(line_order, "N =", counts)
#                      private_alele_plot <-
#                        snps_plot_data %>%
#                        filter(CHROM == current_chr) %>%
#                        ggplot(aes(p_line, POS, color = p_line)) +
#                        geom_point(size = 0.3) + geom_point(size = 0.2, aes(0.5, POS)) +
#                        coord_flip() + labs(x = "Line", y = "Chromossomal Position (Mb)") +
#                        scale_color_discrete(labels = legend, name = "") +
#                        ggtitle(current_chr)
#                      save_plot(paste0("./data/jpegs/snp_", current_chr, ".png"),
#                                private_alele_plot,
#                                base_height = 5, base_aspect_ratio = 2)
#                      return(private_alele_plot)
#                    }, .parallel = TRUE)

u_snp_count <- private_snps %>% count(pu_line, CHROM)
u_snp_count %>% spread(pu_line, n) %>% print(n = 21)

#u_snps_plot_data <-
#  private_snps %>%
#  select(CHROM, POS, pu_line) %>%
#  mutate(POS = POS/1e6)
#psnps_plots = llply(unique(just_snps$CHROM),
#                    function(current_chr)
#                    {
#                      counts = filter(u_snp_count, CHROM == current_chr)
#                      u_line_order = as.character(unique(filter(u_snp_count, CHROM == current_chr)$pu_line))
#                      counts = counts$n[which(u_line_order == filter(u_snp_count,
#                                                                   CHROM == current_chr)$pu_line)]
#                      legend = paste(u_line_order, "N =", counts)
#                      private_alele_plot <-
#                        u_snps_plot_data %>%
#                        filter(CHROM == current_chr) %>%
#                        ggplot(aes(pu_line, POS, color = pu_line)) +
#                        geom_point(size = 0.3) + geom_point(size = 0.2, aes(0.5, POS)) +
#                        coord_flip() + labs(x = "Line", y = "Chromossomal Position (Mb)") +
#                        scale_color_discrete(labels = legend, name = "") +
#                        ggtitle(current_chr)
#                      save_plot(paste0("./data/jpegs/usnps_", current_chr, ".png"),
#                                private_alele_plot,
#                                base_height = 5, base_aspect_ratio = 2)
#                      return(private_alele_plot)
#                    }, .parallel = TRUE)
