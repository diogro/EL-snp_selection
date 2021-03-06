source("./find_all_aleles.R")

all_poly = all_poly %>% arrange(CHROM, POS)

chrom_sizes = all_poly %>%
    group_by(CHROM) %>%
    summarize(start = min(POS), stop = max(POS), diff = (stop - start), min(QUAL) )
chrom_sizes %>%  print(n = 21)

total_size = sum(as.numeric(chrom_sizes$diff[-21]))
total_snps = 50000

snps_per_chrom = floor(as.numeric(chrom_sizes$diff[-21]) * (total_snps/total_size))
names(snps_per_chrom) = chrom_sizes$CHROM[-21]

chunks_per_chrom = floor(snps_per_chrom / 7.3)
names(chunks_per_chrom) = chrom_sizes$CHROM[-21]
sum(chunks_per_chrom)

findChunksStartStop = function(i){
    total = chrom_sizes$diff[i]
    n_chunks = chunks_per_chrom[i]
    chunk_size = floor(total/n_chunks - 1)
    last_chunk = total - (n_chunks - 1) * chunk_size
    chunk_ss = numeric(n_chunks)
    chunk_ss = chunk_ss + chunk_size
    chunk_ss[n_chunks] = last_chunk
    chunk_ss[1] = chunk_ss[1] + chrom_sizes$start[i]
    cumsum(chunk_ss)
}
chunksStartStop = llply(1:20, findChunksStartStop)
names(chunksStartStop) = chrom_sizes$CHROM[-21]

classifyPosition = function(current_chr){
    x = all_poly %>% filter(CHROM == current_chr) %>% select(POS) %>% arrange(POS)
    positions = as.numeric(x$POS)
    css = chunksStartStop[[current_chr]]
    k = 1
    for(i in seq_along(positions)){
        if(positions[i] <= css[k])
            positions[i] = k
        else {k = k+1; positions[i] = k}
    }
    return(data_frame(chunk = positions))
}
chunks = tbl_df(ldply(chrom_sizes$CHROM[-21], classifyPosition, .progress = "text"))
all_poly_noY = filter(all_poly, CHROM != "chrY")
all_poly_noY$CHUNKS = as.numeric(chunks$chunk)
all_poly_list = dlply(all_poly_noY, .(CHROM, CHUNKS), identity)
t(laply(all_poly_list, nrow))

is_good = function(snp, is_private = FALSE, private = c(30, 10), non_private = 20){
    if(!is_private){
        GQ_mask = grepl(".GQ", names(snp))
        all(snp[,GQ_mask] > non_private)
    } else{
        private_mask = grepl(paste0(snp$snp_type, ".GQ"), names(snp))
        GQ_mask = grepl(".GQ", names(snp))
        snp[,private_mask] > private[1] & all(snp[,GQ_mask] > private[2])
    }
}

selectFromAvailable = function(x, available, snps, min_dist, snps_per_class, ...){
    for(current_snp_type in available){
        k = 0
        possible_snps = filter(x, snp_type == current_snp_type) %>% arrange(desc(QUAL))
        possible_snps = bind_rows(filter(possible_snps, aleles != "A/T", aleles != "C/G"),
                                  filter(possible_snps, aleles == "A/T" | aleles == "C/G"))
        for(i in seq(nrow(possible_snps))){
            if(is_good(possible_snps[i,], ...)){
                current_pos = possible_snps[i,]$POS
                dists = abs(current_pos - x$POS)
                dists = dists[dists != 0]
                if(any(dists < min_dist)) next
                new_snp = possible_snps[i,]
                new_snp$priority = ifelse(k == 0 & new_snp$class != "50_50", 1, 0)
                snps = rbind_list(snps, new_snp)
                k = k+1
            }
            if(k > snps_per_class) break
        }
    }
    snps
}

selectSNPchunk = function(x, min_dist = 35, snps_per_class = 4, ...){
    snps = NULL
    available_privates = unique(filter(x, class == "private")$snp_type)
    snps = selectFromAvailable(x, available_privates, snps, min_dist, snps_per_class, is_private = TRUE, ...)
    available_two_four = unique(filter(x, class == "two_four")$snp_type)
    snps = selectFromAvailable(x, available_two_four, snps, min_dist, snps_per_class, ...)
    available_50_50 = unique(filter(x, class == "50_50")$snp_type)
    snps = selectFromAvailable(x, available_50_50, snps, min_dist, snps_per_class, ...)
    available_het = unique(filter(x, class == "private_het")$snp_type)
    snps = selectFromAvailable(x, available_het, snps, min_dist, snps_per_class, ...)
    if(!is.null(snps)){
        snps$nsnps_in_chunk = nrow(snps)
        snps$nsnps_of_type = 0
        snps$qual_rank = 0
        for(type in unique(snps$snp_type)){
            mask = snps$snp_type == type
            n_snps = sum(mask)
            snps$nsnps_of_type[mask] = n_snps
            snps$qual_rank[mask] = 1:n_snps
        }
    }
    snps
}

x = all_poly_list[["chr1.466"]]
current_snp_type = na.omit(unique(x$snp_type))[1]
selectSNPchunk(x) %>% select(POS, QUAL, TYPE, aleles, class, snp_type, priority, nsnps_in_chunk, nsnps_of_type, qual_rank)

new_all_poly_noY = (all_poly_noY)
all_poly_list = dlply(new_all_poly_noY, .(CHROM, CHUNKS), identity)
selected_snps_noY  = tbl_df(ldply(all_poly_list, selectSNPchunk, .parallel = TRUE)) %>% arrange(CHROM, POS)
selected_snps_noY %>% count(aleles)

selected_dist = diff(selected_snps_noY$POS)
quantile(selected_dist[selected_dist>0], seq(0, 1, 0.01))
too_close = which(selected_dist == 11 | selected_dist == 20)
selected_dist[too_close]
selected_snps = selected_snps_noY[-c(too_close, too_close+1),]

chromY = all_poly %>% filter(CHROM == "chrY")
snp_y_mask = names(chromY)
snp_y_mask = (grepl("_1", snp_y_mask) | grepl("_2", snp_y_mask))
chromY[,snp_y_mask] %>% print(n=100)
chromY %>% count(TYPE, class)
y_snps = selectSNPchunk(chromY, min_dist = 1, 100, private = c(0, 0), non_private = 0) %>% filter(class != "private_het") %>% arrange(POS)
y_snps %>% select(POS, QUAL, TYPE, class, snp_type, priority) %>% print(n = 50)
selected_snps =bind_rows(selected_snps, y_snps)

sel_usnp_count <- selected_snps %>% count(snp_type, CHROM)
sel_usnp_count %>% spread(snp_type, n) %>% print(n = 21)
sel_usnp_quality <- selected_snps %>% group_by(CHROM) %>% summarize(min(QUAL), quantile(QUAL, 0.001), median(QUAL))
print(sel_usnp_quality, n = 21)

u_snps_plot_data <-
    selected_snps %>%
    select(CHROM, POS, snp_type) %>%
    mutate(POS = POS/1e6)

#snp_order = as.character(unique(sel_usnp_count$snp_type))
#snp_order = snp_order[order(nchar(snp_order))]
#u_snps_plot_data$snp_type = factor(u_snps_plot_data$snp_type, levels = snp_order)
#psnps_plots = llply(unique(selected_snps$CHROM),
                    #function(current_chr)
                    #{
                        #counts = filter(sel_usnp_count, CHROM == current_chr)
                        #u_line_order = as.character(unique(filter(sel_usnp_count, CHROM == current_chr)$snp_type))
                        #u_line_order = u_line_order[order(nchar(u_line_order))]
                        #counts = counts$n[match(u_line_order, counts$snp_type)]
                        #legend = paste(u_line_order, "N =", counts)
                        #private_alele_plot <-
                            #u_snps_plot_data %>%
                            #filter(CHROM == current_chr) %>%
                            #ggplot(aes(snp_type, POS, color = snp_type)) +
                            #geom_point(size = 0.2) + geom_point(size = 0.2, aes(0.5, POS)) +
                            #coord_flip() + labs(x = "Line", y = "Chromossomal Position (Mb)") +
                            #scale_color_discrete(labels = legend, name = "") +
                            #ggtitle(current_chr)
                        #save_plot(paste0("./data/jpegs/selected_usnps_", current_chr, ".png"),
                                  #private_alele_plot,
                                  #base_height = 6.5, base_aspect_ratio = 2.3)
                        #return(private_alele_plot)
                    #}, .parallel = TRUE)

selected_table = selected_snps %>% 
    unite("A13", A13_1, A13_2, sep = "/") %>%
    unite("A31", A31_1, A31_2, sep = "/") %>%
    unite("A41", A41_1, A41_2, sep = "/") %>%
    unite("A23", A23_1, A23_2, sep = "/") %>%
    unite("A22", A22_1, A22_2, sep = "/") %>%
    unite("A42", A42_1, A42_2, sep = "/") %>%
    select(-is_bialelic, -TYPE)
write_csv(selected_table, "./data/selected_snps.csv")
