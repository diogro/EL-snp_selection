source("./find_private_aleles.R")

chrom_sizes = all_snps %>%
    group_by(CHROM) %>%
    summarize(start = min(POS), stop = max(POS), diff = (stop - start), min(QUAL) )
chrom_sizes %>%  print(n = 21)

total_size = sum(as.numeric(chrom_sizes$diff[-21]))
total_snps = 50000

snps_per_chrom = floor(as.numeric(chrom_sizes$diff[-21]) * (total_snps/total_size))
names(snps_per_chrom) = chrom_sizes$CHROM[-21]

chunks_per_chrom = floor(snps_per_chrom / 12)
names(chunks_per_chrom) = chrom_sizes$CHROM[-21]

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
    x = all_snps %>% filter(CHROM == current_chr) %>% select(POS)
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
all_snps_noY = filter(all_snps, CHROM != "chrY")
all_snps_noY$CHUNKS = as.numeric(chunks$chunk)
all_snps_list = dlply(all_snps_noY, .(CHROM, CHUNKS), identity)
t(laply(all_snps_list, nrow))

mean_qual = ddply(all_snps_noY[,c("CHROM", "CHUNKS", names(all_snps_list[[1]])[grep(".GQ", names(all_snps_list[[1]]))])],
      .(CHROM, CHUNKS), function(x) colMeans(x[,-c(1, 2)]))
min_qual = ddply(mean_qual, .(CHROM, CHUNKS), function(x) min(x[,-c(1, 2)]))
qual_plot = ggplot(min_qual, aes(CHUNKS, V1, group = CHROM, color= CHROM)) + geom_point()
save_plot(paste0("./data/jpegs/qual_chunk.png"),
          qual_plot,
          base_height = 5, base_aspect_ratio = 2)

is_good_p = function(snp){
    private_mask = grepl(paste0(snp$p_line, ".GQ"), names(snp))
    GQ_mask = grepl(".GQ", names(snp))
    snp[,private_mask] > 30 & all(snp[,GQ_mask] > 10)
}
is_good = function(snp){
    GQ_mask = grepl(".GQ", names(snp))
    all(snp[,GQ_mask] > 20)
}
selectSNPchunk = function(x){
    available_privates = na.omit(unique(x$p_line))
    psnp_df = rbind_all(llply(available_privates, function(current_line) {
                        possible_snps = filter(x, p_line == current_line) %>% arrange(desc(QUAL))
                        for(i in seq(nrow(possible_snps))){
                            snp = possible_snps[i,]
                            if(is_good_p(snp)){ break
                            } else snp = NULL
                        }
                        snp
                        }
                  )
             )
    available_usefull = na.omit(unique(x$pu_line[!x$is_private]))
    usnp_df = rbind_all(llply(available_usefull, function(current_line) {
                        possible_snps = filter(x, pu_line == current_line) %>% arrange(desc(QUAL))
                        for(i in seq(nrow(possible_snps))){
                            snp = possible_snps[i,]
                            if(is_good(snp)){ break
                            } else snp = NULL
                        }
                        snp
                        }
                  )
             )
    rbind_list(psnp_df, usnp_df)
}
selectSNPchunk(x)
new_all_snps_noY = (all_snps_noY)
all_snps_list = dlply(new_all_snps_noY, .(CHROM, CHUNKS), identity)
selected_snps = tbl_df(ldply(all_snps_list, selectSNPchunk, .parallel = TRUE))

new_all_snps_noY = setdiff(new_all_snps_noY, selected_snps)
all_snps_list = dlply(new_all_snps_noY, .(CHROM, CHUNKS), identity)
selected_snps = bind_rows(selected_snps, tbl_df(ldply(all_snps_list, selectSNPchunk, .parallel = TRUE)))

sel_usnp_count <- selected_snps %>% count(pu_line, CHROM)
sel_usnp_count %>% spread(pu_line, n) %>% print(n = 21)

u_snps_plot_data <-
  selected_snps %>%
  select(CHROM, POS, pu_line) %>%
  mutate(POS = POS/1e6)
psnps_plots = llply(unique(selected_snps$CHROM),
                    function(current_chr)
                    {
                      counts = filter(sel_usnp_count, CHROM == current_chr)
                      u_line_order = as.character(unique(filter(sel_usnp_count, CHROM == current_chr)$pu_line))
                      counts = counts$n[which(u_line_order == filter(sel_usnp_count,
                                                                   CHROM == current_chr)$pu_line)]
                      legend = paste(u_line_order, "N =", counts)
                      private_alele_plot <-
                        u_snps_plot_data %>%
                        filter(CHROM == current_chr) %>%
                        ggplot(aes(pu_line, POS, color = pu_line)) +
                        geom_point(size = 0.3) + geom_point(size = 0.2, aes(0.5, POS)) +
                        coord_flip() + labs(x = "Line", y = "Chromossomal Position (Mb)") +
                        scale_color_discrete(labels = legend, name = "") +
                        ggtitle(current_chr)
                      save_plot(paste0("./data/jpegs/selected_usnps_", current_chr, ".png"),
                                private_alele_plot,
                                base_height = 5, base_aspect_ratio = 2)
                      return(private_alele_plot)
                    }, .parallel = TRUE)
