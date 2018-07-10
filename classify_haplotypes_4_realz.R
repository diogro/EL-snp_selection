source("read_phenotypes.R")

if(!require(vcfR)){install.packages("vcfR"); library(vcfR)}
if(!require(stringr)){install.packages("stringr"); library(stringr)}
vcf <- read.vcfR("./data/imputed.vcf.vcf.gz")
raw_gen = tbl_df(cbind(vcf@fix, vcf@gt))

full_id = colnames(vcf@gt)[-1]
IDs = ldply(full_id, function(x) c(x, unlist(strsplit(x, "_"))))
colnames(IDs) = c("full", "ID")
rm(vcf)

gen = dplyr::select(raw_gen, CHROM, POS, ID, everything()) %>%
    dplyr::select(-REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT) %>%
    rename(chr = CHROM,
           pos = POS) %>%
    arrange(chr, pos) %>%
    filter(chr != 21)
gen$chr = as.numeric(gen$chr)
gen$pos = as.numeric(gen$pos)
rm(raw_gen)

## From http://churchill-lab.jax.org/mousemapconverter
## NCBI Build 37 bp -> Sex-Averaged cM Cox
positions = read.table("./data/marker_positions.txt", header = FALSE, stringsAsFactors = FALSE) %>% tbl_df
names(positions) = c("chr2", "pos", "chr", "gpos")
positions$chr = as.integer(positions$chr)
tail(positions$chr)
gen = inner_join(positions, gen, by = c("chr", "pos")) %>%
  dplyr::select(ID, chr, everything(), -chr2)
rownames(gen) = gen$ID
chr_maps = dlply(gen, .(chr), function(x) x[,1:4])

full_snped = inner_join(full_data, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Foster_litter_size_weaning, Final_weight, contains("growth"))

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Foster_litter_size_weaning, Final_weight, contains("growth"))

f5_snped = inner_join(full_data_F5, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Foster_litter_size_weaning, Final_weight, contains("growth"))

f1_snped = inner_join(full_data_F1, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Foster_litter_size_weaning, Final_weight, contains("growth"))

pat_snped = inner_join(full_data_Strain, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Foster_litter_size_weaning, Final_weight, contains("growth"))

.n = function(x) as.numeric(factor(x, levels = c("M", "F")))
pedigree = as.data.frame(read.csv("./data/Intercross_pedigree2.csv")) %>% dplyr::rename(id = animal) %>% orderPed
full_data$ID = as.numeric(full_data$ID)
ped2 = (left_join(dplyr::rename(pedigree, ID = id), dplyr::select(full_data, ID, Sex), by = "ID"))
full_data$ID = as.character(full_data$ID)
missing = full_snped[!full_snped$ID %in% ped2$ID,c("ID", "Mat_ID", "Pat_ID", "Sex")]
names(missing) = names(ped2)
ped2 = rbind(ped2, missing) %>% orderPed
missing_sire = data.frame(ID = unique(na.omit(ped2$sire[!ped2$sire %in% ped2$ID])), dam = NA, sire = NA, Sex = "M")
ped2 = rbind(ped2, missing_sire) %>% orderPed
missing_dam = data.frame(ID = unique(na.omit(ped2$dam[!ped2$dam %in% ped2$ID])), dam = NA, sire = NA, Sex = "F")
ped2 = rbind(ped2, missing_dam) %>% orderPed
for(i in 1:nrow(ped2)){
  if(is.na(ped2$Sex[i])){
    if(ped2$ID[i] %in% ped2$dam) ped2$Sex[i] = "F"
    else if(ped2$ID[i] %in% ped2$sire) ped2$Sex[i] = "M"
    else ped2$Sex[i] = "M"
  }
}

A = 2*kinship(pedigree)
A = A + diag(nrow(A)) * 1e-4
ids = as.character(f6_snped$ID)
Af6 = (A[ids,ids])
colnames(Af6) = rownames(Af6) = f6_snped$ID
load("./data/gemma_relatedness.Rdata")
#pedAll <- pedigree(id=ped2$ID,
#dadid=ped2$sire, momid=ped2$dam,
#sex=ped2$Sex)

getFounders = function(current_line){
  gen_cols = IDs$full[match(as.character(dplyr::filter(pat_snped, Strain == current_line)$ID), IDs$ID)]
  current_line_gen = gen[,gen_cols]
  current_line_gen = t(laply(current_line_gen, substring, 1, 3))
  colnames(current_line_gen) = paste(current_line, 1:length(gen_cols), sep = "_")
  rownames(current_line_gen) = gen$ID
  current_line_gen
}
createFounderDict = function(i){
    f_dict = vector("list", 2)
    f_dict[[1]] = names(founders)[sapply(founders, function(x) any(str_detect(x[i,], "1")))]
    f_dict[[2]] = names(founders)[sapply(founders, function(x) any(str_detect(x[i,], "0")))]
    names(f_dict) = c("1", "0")
    return(f_dict)
}
classifyInd = function(ID, join = TRUE){
  print(ID)
  out = vector("list", 2)
  gen = as.data.frame(gen)
  tape = c(1, 3)
  for(k in 1:2){
    single_ind = Map(`[[`, founders_dict, substring(gen[,ID], tape[k],tape[k]))
    single_ind_chr = vector("list", 20)
    for(i in 1:20) single_ind_chr[[i]] = single_ind[gen$chr == i]
    haplotypes = tibble(chr = numeric(),
                        seg = numeric(),
                        start = numeric(),
                        finish = numeric(),
                        len = numeric(),
                        strain = character(),
                        nstrain = numeric())
    num_seg = 0
    for(chr in 1:20){
      seg_start = 1
      num_seg = num_seg + 1
      seg = 1
      plausible_lines = line_order
      n_marker = length(single_ind_chr[[chr]])
      for(i in seq_along(single_ind_chr[[chr]])){
        plausible_lines_new = intersect(plausible_lines, single_ind_chr[[chr]][[i]])
        if(length(plausible_lines_new) == 0 & i != n_marker){
          seg_start = i
          seg = seg + 1
          num_seg = num_seg + 1
          plausible_lines = single_ind_chr[[chr]][[i]]
          line = c(chr,
                   "R",
                   seg_start, i, 1,
                   paste(plausible_lines, collapse = ":"), length(plausible_lines))
          #print(line)
          haplotypes[num_seg,] = line
          num_seg = num_seg + 1
          plausible_lines = line_order
        } else{
          if(length(plausible_lines_new) > 0) plausible_lines = plausible_lines_new
          line = c(chr,
                   seg,
                   seg_start, i, i - seg_start,
                   paste(plausible_lines, collapse = ":"), length(plausible_lines))
          #print(line)
          haplotypes[num_seg,] = line
        }
      }
    }
    out[[k]] = haplotypes
  }
  if(join){
    out = lapply(out, fixRecombinations)
    out = lapply(out, addPositions)
  }
  out
}
fixRecombinations = function(x){
  while(TRUE){
    recon = which(x$seg == "R")[1]
    #if(is.na(recon)) recon = which(x$len == 1)[1]
    if(is.na(recon)) break

    ## Check if previous and next segments can be joined
    chr = as.numeric(x$chr[recon])
    start = as.numeric(x$start[recon-1])
    middle = as.numeric(x$start[recon])
    end = as.numeric(x$finish[recon+1])
    consensus = intersect(unlist(strsplit(x$strain[recon-1], ":")),
                          unlist(strsplit(x$strain[recon+1], ":")))
    if(length(consensus) > 0){
      seg = x$seg[recon-1]
      x = x[-c(recon-1, recon, recon+1),]
      line = c(chr,
               seg,
               start, end, end - start,
               paste(consensus, collapse = ":"), length(consensus))
      x[nrow(x)+1,] = line
      x = arrange(x, as.numeric(chr), as.numeric(start))
      ## Check if next segment can be joined
    } else{
      consensus = intersect(unlist(strsplit(x$strain[recon], ":")),
                            unlist(strsplit(x$strain[recon+1], ":")))
      if(length(consensus) > 0){
        seg = x$seg[recon+1]
        x = x[-c(recon, recon+1),]
        line = c(chr,
                 seg,
                 middle, end, end - middle,
                 paste(consensus, collapse = ":"), length(consensus))
        x[nrow(x)+1,] = line
        x = arrange(x, as.numeric(chr), as.numeric(start))
      } else{
        seg = x$seg[recon+1]
        next_seg_strain = x$strain[recon+1]
        x = x[-c(recon, recon+1),]
        line = c(chr,
                 seg,
                 middle, end, end - middle,
                 next_seg_strain, length(next_seg_strain))
        x[nrow(x)+1,] = line
        x = arrange(x, as.numeric(chr), as.numeric(start))
      }
    }
  }
  x
}
addPositions = function(x){
  x$start_pos = unlist(sapply(1:nrow(x), function(i) chr_maps[[x$chr[i]]][x$start[i],"pos"]))
  x$finish_pos = unlist(sapply(1:nrow(x), function(i) chr_maps[[x$chr[i]]][x$finish[i],"pos"]))
  x$len_bp = x$finish_pos - x$start_pos
  x
}
plotHaplotype = function(hap){
  ids = c(paste(1, hap[[1]]$chr, hap[[1]]$seg, sep = "."),
          paste(2, hap[[2]]$chr, hap[[2]]$seg, sep = "."))

  values <- data.frame(
    id = ids,
    value = c(hap[[1]]$strain, hap[[2]]$strain)
  )
  ycoord = c()
  xcoord = c()
  for(i in 1:nrow(hap[[1]])) {
    ycoord = c(ycoord, hap[[1]]$start_pos[i], hap[[1]]$finish_pos[i], hap[[1]]$finish_pos[i], hap[[1]]$start_pos[i],hap[[1]]$start_pos[i])
    xcoord = c(xcoord,
               as.numeric(hap[[1]]$chr)[i] - 0.4,
               as.numeric(hap[[1]]$chr)[i] - 0.4,
               as.numeric(hap[[1]]$chr)[i] - 0.05,
               as.numeric(hap[[1]]$chr)[i] - 0.05,
               as.numeric(hap[[1]]$chr)[i] - 0.4)
  }
  for(i in 1:nrow(hap[[2]])) {
    ycoord = c(ycoord, hap[[2]]$start_pos[i], hap[[2]]$finish_pos[i], hap[[2]]$finish_pos[i], hap[[2]]$start_pos[i],hap[[2]]$start_pos[i])
    xcoord = c(xcoord,
               as.numeric(hap[[2]]$chr)[i] + 0.05,
               as.numeric(hap[[2]]$chr)[i] + 0.05,
               as.numeric(hap[[2]]$chr)[i] + 0.4,
               as.numeric(hap[[2]]$chr)[i] + 0.4,
               as.numeric(hap[[2]]$chr)[i] + 0.05)
  }
  bar_positions <- data.frame(
    id = rep(ids, each = 5),
    x = xcoord,
    y = as.numeric(ycoord)
  )
  datapoly <- merge(values, bar_positions, by = c("id"))
  p = ggplot(datapoly, aes(x = x, y = y)) +
    geom_polygon(aes(fill = value, group = id), color = "black") +
    scale_x_continuous(breaks = 1:20, labels = c(1:19, "X"))
  return(p)
}

founders = llply(line_order, getFounders)
names(founders) = line_order

founders_dict = lapply(1:nrow(founders[[1]]), createFounderDict)
names(founders_dict) = gen$ID

# registerDoMC(4)
# f6_haplotypes = llply(f6_snped$ID, classifyInd, .parallel = TRUE)
# names(f6_haplotypes) = f6_snped$ID
# f5_haplotypes = llply(f5_snped$ID, classifyInd, .parallel = TRUE)
# names(f5_haplotypes) = f5_snped$ID
# f1_haplotypes = llply(f1_snped$ID, classifyInd, .parallel = TRUE)
# names(f1_haplotypes) = f1_snped$ID
# save(founders, founders_dict, f6_haplotypes, f5_haplotypes, f1_haplotypes, file = "./data/infered_haplotypes.Rdata")
load("./data/infered_haplotypes.Rdata")

# for(i in seq_along(f1_snped$ID)){
#  p = plotHaplotype(f1_haplotypes[[i]]) + labs(x = "chromossome", y = "Position (bp)") + ggtitle(paste("F1",f1_snped$ID[[i]], f1_snped$Sex[i]))
#  save_plot(paste0("./data/haplotype_inference_figures/f1_haplotype_ID-", f1_snped$ID[i], ".png"), p, base_height = 5, base_aspect_ratio = 2)
# }
# for(i in seq_along(f6_snped$ID)){
#   p = plotHaplotype(f6_haplotypes[[i]]) + labs(x = "chromossome", y = "Position (bp)") + ggtitle(paste("F6",f6_snped$ID[[i]], f6_snped$Sex[i]))
#   save_plot(paste0("./data/haplotype_inference_figures/f6_haplotype_ID-", f6_snped$ID[i], ".png"), p, base_height = 5, base_aspect_ratio = 2)
# }
# for(i in seq_along(f5_snped$ID)){
#   p = plotHaplotype(f5_haplotypes[[i]]) + labs(x = "chromossome", y = "Position (bp)") + ggtitle(paste("F5",f5_snped$ID[[i]], f5_snped$Sex[i]))
#   save_plot(paste0("./data/haplotype_inference_figures/f5_haplotype_ID-", f5_snped$ID[i], ".png"), p, base_height = 5, base_aspect_ratio = 2)
# }

# compareChr = function(x, y) sum(out <- sapply(Map(intersect, lapply(x, strsplit, ":"), lapply(y, strsplit, ":")), length))/length(out)
#
# ID = "3486"
# sire_ID = as.character(pedigree[pedigree$id == as.numeric(ID),"sire"])
# dam_ID = as.character(pedigree[pedigree$id == as.numeric(ID),"dam"])
# i = 2
# for(i in 1:20){
#   sire_hap_1 = filter(f5_haplotypes[[sire_ID]][[1]], chr == i)
#   sire_hap_2 = f5_haplotypes[[sire_ID]][[2]] %>% filter(chr == i)
#   dam_hap_1 = f5_haplotypes[[dam_ID]][[1]] %>% filter(chr == i)
#   dam_hap_2 = f5_haplotypes[[dam_ID]][[2]] %>% filter(chr == i)
#   offspring_hap_1 = filter(f6_haplotypes[[ID]][[1]], chr == i)
#   offspring_hap_2 = filter(f6_haplotypes[[ID]][[2]], chr == i)
#
#   sire_seq_1 = rep(sire_hap_1$strain, as.numeric(sire_hap_1$len) + 1)
#   sire_seq_2 = rep(sire_hap_2$strain, as.numeric(sire_hap_2$len) + 1)
#   dam_seq_1 = rep(dam_hap_1$strain, as.numeric(dam_hap_1$len) + 1)
#   dam_seq_2 = rep(dam_hap_2$strain, as.numeric(dam_hap_2$len) + 1)
#   offspring_seq_1 = rep(offspring_hap_1$strain, as.numeric(offspring_hap_1$len) + 1)
#   offspring_seq_2 = rep(offspring_hap_2$strain, as.numeric(offspring_hap_2$len) + 1)
#
#   sapply(list(sire_seq_1, sire_seq_2, dam_seq_1, dam_seq_2), compareChr, offspring_seq_1)
#   sapply(list(sire_seq_1, sire_seq_2, dam_seq_1, dam_seq_2), compareChr, offspring_seq_2)
# }
# sire_hap_1 %>% print(n = nrow(.))

createModelMatrix = function(x){
  strain_seq_1 = rep(x[[1]]$strain, as.numeric(x[[1]]$len) + 1)
  strain_seq_2 = rep(x[[2]]$strain, as.numeric(x[[2]]$len) + 1)
  singleMarkerLine = function(x){
    marker_strains = strsplit(x, ":")[[1]]
    as.numeric(line_order %in% marker_strains)/length(marker_strains)
  }
  model = laply(strain_seq_1, singleMarkerLine) + laply(strain_seq_2, singleMarkerLine)
  colnames(model) = line_order
  rownames(model) = gen$ID
  model
}
# f6_model_matrices = laply(f6_haplotypes, createModelMatrix, .parallel = TRUE)
# f6_model_matrices = aperm(f6_model_matrices, c(1, 3, 2))
# dimnames(f6_model_matrices)[[1]] = f6_snped$ID
# save(f6_model_matrices, file = "./data/f6_model_matrices.Rdata")
load("./data/f6_model_matrices.Rdata")
