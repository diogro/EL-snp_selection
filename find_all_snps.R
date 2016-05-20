source("./read_raw_snp_data.R")

get_private = function(x, tx = table(x)){
  if(any(tx == 2)){
    rare = names(tx)[which(min(tx) == tx)]
    lines = names(x)[x == rare]
    line = unique(laply(strsplit(lines, "_"), `[`, 1))
    if(length(line) == 1) return(line)
    else return(NULL)
  }
  else return(NULL)
}

get_two_four = function(x, tx = table(x)){
  if(any(tx == 4)){
    rare = names(tx)[which(min(tx) == tx)]
    lines = names(x)[x == rare]
    line = unique(laply(strsplit(lines, "_"), `[`, 1))
    line = line_order[line_order %in% line]
    if(length(line) == 2) return(paste(line, collapse = ""))
    else return(NULL)
  }
  else return(NULL)
}

get_50_50 = function(x, tx = table(x)){
    if(all(tx == 6)){
        alele_1 = names(tx)[1]
        alele_2 = names(tx)[2]
        lines_1 = names(x)[x == alele_1]
        lines_2 = names(x)[x == alele_2]
        lines_1 = unique(laply(strsplit(lines_1, "_"), `[`, 1))
        lines_2 = unique(laply(strsplit(lines_2, "_"), `[`, 1))
        if(length(lines_1) == 3 & length(lines_2) == 3){
            if("A13" %in% lines_1){
                lines_1 = line_order[line_order %in% lines_1]
                return(paste(lines_1, collapse = ""))
            } else{
                lines_2 = line_order[line_order %in% lines_2]
                return(paste(lines_2, collapse = ""))
            }
        }
        else return(NULL)
    }
    else return(NULL)
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
    snp_kind = NULL
    if(!is.null(snp_kind <- get_private(x, tx))) return(c("private", snp_kind))
    if(!is.null(snp_kind <- get_two_four(x, tx))) return(c("two_four", snp_kind))
    if(!is.null(snp_kind <- get_50_50(x, tx))) return(c("50_50", snp_kind))
    if(any(tx == 1)){
        if(get_A22_p_het(x)) return(c("private_het", "A22_pri_het"))
    } else if(any(tx == 3)){
        if(get_A22_het_p_42(x)) return(c("private_het", "A22_het_p_A42"))
        else if(get_A22_het_p_23(x)) return(c("private_het", "A22_het_p_A23"))
    }
    return(c(NA, NA))
}

#is_snp = tbl_df(ldply(snp_array_list,
                      #function(snp_array) adply(snp_array, 1, is_any_kind),
                      #.parallel = TRUE))
#save(is_snp, file = "./data/all_snp_position.Rdata")
load("./data/all_snp_position.Rdata")

line_order = c("A13", "A31", "A41", "A23", "A22", "A42")[6:1]
bialelic_snps = mutate(bialelic_snps,
                       class = is_snp$V1,
                       snp_type = is_snp$V2) %>% select(-CHUNKS)

not_bialelic_snps = all_poly %>%
    filter(!is_bialelic) %>%
    mutate(class = NA,
           snp_type = NA)

all_poly = bind_rows(bialelic_snps, not_bialelic_snps) %>% arrange(CHROM, POS)

private_snps = all_poly %>% filter(class == "private")

p_snp_count <-
  all_poly %>%
  filter(class == "private") %>%
  count(snp_type, CHROM)
(p_snp_count) %>% spread(snp_type, n) %>% print(n = 21)

#snps_plot_data <-
  #all_poly %>%
  #filter(class == "private") %>%
  #select(CHROM, POS, snp_type) %>%
  #mutate(POS = POS/1e6)
#psnps_plots = llply(unique(just_snps$CHROM),
                    #function(current_chr)
                    #{
                      #counts = filter(p_snp_count, CHROM == current_chr)
                      #counts = counts$n[which(line_order == filter(p_snp_count,
                                                                   #CHROM == current_chr)$snp_type)]
                      #legend = paste(line_order, "N =", counts)
                      #private_alele_plot <-
                        #snps_plot_data %>%
                        #filter(CHROM == current_chr) %>%
                        #ggplot(aes(snp_type, POS, color = snp_type)) +
                        #geom_point(size = 0.3) + geom_point(size = 0.2, aes(0.5, POS)) +
                        #coord_flip() + labs(x = "Line", y = "Chromossomal Position (Mb)") +
                        #scale_color_discrete(labels = legend, name = "") +
                        #ggtitle(current_chr)
                      #save_plot(paste0("./data/jpegs/snp_", current_chr, ".png"),
                                #private_alele_plot,
                                #base_height = 5, base_aspect_ratio = 2)
                      #return(private_alele_plot)
                    #}, .parallel = TRUE)

u_snp_count <- all_poly %>% filter(!is.na(class)) %>% count(snp_type, CHROM)
u_snp_count %>% spread(snp_type, n) %>% print(n = 21)

#u_snps_plot_data <-
  #all_poly %>% filter(!is.na(class)) %>%
  #select(CHROM, POS, snp_type) %>%
  #mutate(POS = POS/1e6)
#snp_order = as.character(unique(u_snp_count$snp_type))
#snp_order = snp_order[order(nchar(snp_order))]
#u_snps_plot_data$snp_type = factor(u_snps_plot_data$snp_type, levels = snp_order)
#current_chr = "chr16"
#psnps_plots = llply(unique(just_snps$CHROM),
                    #function(current_chr)
                    #{
                      #counts = filter(u_snp_count, CHROM == current_chr)
                      #u_line_order = as.character(unique(filter(u_snp_count, CHROM == current_chr)$snp_type))
                      #u_line_order = u_line_order[order(nchar(u_line_order))]
                      #counts = counts$n[match(u_line_order, counts$snp_type)]
                      #legend = paste(u_line_order, "n =", counts)
                      #private_alele_plot <-
                        #u_snps_plot_data %>% filter(CHROM == current_chr) %>%
                        #ggplot(aes(snp_type, POS, color = snp_type)) +
                        #geom_point(size = 0.1) + geom_point(size = 0.1, aes(0.5, POS)) +
                        #coord_flip() + labs(x = "Type", y = "Chromossomal position (mb)") +
                        #scale_color_discrete(labels = legend, name = "") +
                        #ggtitle(current_chr)
                      #save_plot(paste0("./data/jpegs/usnps_", current_chr, ".png"),
                                #private_alele_plot,
                                #base_height = 6, base_aspect_ratio = 2.2)
                      #return(private_alele_plot)
                    #}, .parallel = TRUE)
