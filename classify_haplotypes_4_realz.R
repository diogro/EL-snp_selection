source("read_phenotypes.R")

if(!require(rcfR)){install.packages("vcfR"); library(vcfR)}
if(!require(stringr)){install.packages("stringr"); library(stringr)}
vcf <- read.vcfR("./data/imputed.vcf.vcf.gz")
raw_gen = tbl_df(cbind(vcf@fix, vcf@gt))

full_id = colnames(vcf@gt)[-1]
IDs = ldply(full_id, function(x) c(x, unlist(strsplit(x, "_"))))
colnames(IDs) = c("full", "ID")

gen = dplyr::select(raw_gen, CHROM, POS, ID, everything()) %>%
    dplyr::select(-REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT) %>%
    rename(chr = CHROM,
           pos = POS) %>%
    arrange(chr, pos) %>%
    filter(chr != 21)
gen$chr = as.numeric(gen$chr)
gen$pos = as.numeric(gen$pos)

## From http://churchill-lab.jax.org/mousemapconverter
## NCBI Build 37 bp -> Sex-Averaged cM Cox
positions = read.table("./data/marker_positions.txt", header = FALSE, stringsAsFactors = FALSE) %>% tbl_df
names(positions) = c("chr2", "pos", "chr", "gpos")
positions$chr = as.integer(positions$chr)
tail(positions$chr)
gen = inner_join(positions, gen, by = c("chr", "pos")) %>%
  dplyr::select(ID, chr, everything(), -chr2)
rownames(gen) = gen$ID
full_snped = inner_join(full_data, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID, Final_weight)


p = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID)

f5_snped = inner_join(full_data_F5, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID)

f1_snped = inner_join(full_data_F1, IDs, by = "ID") %>%
  select(ID, Strain, Litter_ID_new:Mat_ID)

pat_snped = inner_join(full_data_Strain, IDs, by = "ID") %>%
  select(ID, Strain, Litter_ID_new:Mat_ID)

nrow(pat_snped)

current_line = line_order[5]
getFounders = function(current_line){
  gen_cols = IDs$full[match(as.character(dplyr::filter(pat_snped, Strain == current_line)$ID), IDs$ID)]
  current_line_gen = gen[,gen_cols]
  current_line_gen = t(laply(current_line_gen, substring, 1, 3))
  colnames(current_line_gen) = paste(current_line, 1:length(gen_cols), sep = "_")
  rownames(current_line_gen) = gen$ID
  current_line_gen
}
founders = llply(line_order, getFounders)
names(founders) = line_order
createFounderDict = function(i){
    f_dict = vector("list", 2)
    f_dict[[1]] = names(founders)[sapply(founders, function(x) any(str_detect(x[i,], "1")))]
    f_dict[[2]] = names(founders)[sapply(founders, function(x) any(str_detect(x[i,], "0")))]
    names(f_dict) = c("1", "0")
    return(f_dict)
}
founders_dict = lapply(1:nrow(founders[[1]]), createFounderDict)
names(founders_dict) = gen$ID

ID = "3033"
classifyInd = function(ID){
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
      for(i in seq_along(single_ind_chr[[chr]])){
        plausible_lines = intersect(plausible_lines, single_ind_chr[[chr]][[i]])
        if(length(plausible_lines) == 0){
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
  lapply(out, fixRecombinations)
}
x = out[[2]]
fixRecombinations = function(x){
  while(TRUE){
    recon = which(x$seg == "R")[1]
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
        x = x[-c(recon, recon+1),]
        line = c(chr,
                 seg,
                 middle, end, end - middle,
                 x$strain[recon+1], length(x$strain[recon+1]))
        x[nrow(x)+1,] = line
        x = arrange(x, as.numeric(chr), as.numeric(start))
      }
    }
  }
  x
}
print(x, n= nrow(x))

hap = classifyInd("3034")
print(hap[[1]], n = nrow(hap[[1]]))
print(hap[[2]], n = nrow(hap[[2]]))
f6_haplotypes = llply(f6_snped$ID, classifyInd, .progress = "text")
