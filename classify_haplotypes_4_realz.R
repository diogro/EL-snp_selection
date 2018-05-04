source("read_phenotypes.R")

if(!require(rcfR)){install.packages("vcfR"); library(vcfR)}
if(!require(stringr)){install.packages("stringr"); library(stringr)}
vcf <- read.vcfR("./data/plink_files/imputed.vcf")
raw_gen = tbl_df(cbind(vcf@fix, vcf@gt))
  
full_id = colnames(vcf@gt)[-1]
IDs = ldply(full_id, function(x) c(x, unlist(strsplit(x, "_"))))
colnames(IDs) = c("full", "fam", "ID")

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

full_snped = inner_join(full_data, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID, Final_weight)

f6_snped = inner_join(full_data_F6, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID)

f5_snped = inner_join(full_data_F5, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID)

f1_snped = inner_join(full_data_F1, IDs, by = "ID") %>%
  select(ID, Litter_ID_new:Mat_ID)

pat_snped = inner_join(full_data_Strain, IDs, by = "ID") %>%
  select(ID, Strain, Litter_ID_new:Mat_ID)

nrow(pat_snped)

current_line = line_order[5]
getFounders = function(current_line){
  gen_cols = IDs$full[match(as.character(dplyr::filter(pat_snped, Strain == current_line)$ID), IDs$ID)]
  current_line_gen = gen %>% dplyr::select(ID, chr, pos, gpos, gen_cols)
  qtl2convert::find_consensus_geno(as.matrix(current_line_gen[,-c(1:4)]))
}
founders = gsub("/", "", laply(line_order, getFounders))
dimnames(founders)[[1]] = line_order
dimnames(founders)[[2]] = gen$ID
founders = t(founders)
founders = substring(founders, 1, 3)

nrow(founders[apply(is.na(founders), 1, any),])
colnames(founders)

x = founders[1,]
createFounderDict = function(x){
    f_dict = vector("list", 2)
    f_dict[[1]] = names(x)[str_detect(x, "1")]
    f_dict[[2]] = names(x)[str_detect(x, "0")]
    names(f_dict) = c("1", "0")
    return(f_dict)
}
founders_dict = apply(founders, 1, createFounderDict)

single_ind = Map(`[[`, founders_dict, substring(gen$`2_127`, 1, 1))

table(laply(single_ind, length))

haplotypes = tibble(start = numeric(), finish = numeric(), len = numeric(), strain = character(), nstrain = numeric())
k = 1
l = 0
plausible_lines = line_order
for(i in seq_along(single_ind)){
    plausible_lines_new = plausible_lines[which(plausible_lines %in% single_ind[[i]])]
    if(length(plausible_lines_new) != 0){
        plausible_lines = plausible_lines_new
    } else {
        l = l + 1
        if(l == 10){
            haplotypes[nrow(haplotypes)+1,] = c(k, i, i-k, paste(plausible_lines, collapse = ":"), length(plausible_lines))
            k = i
            l = 0
            plausible_lines = line_order
            plausible_lines_new = plausible_lines[which(plausible_lines %in% single_ind[[i]])]
        }
    }
}
print(haplotypes, n = nrow(haplotypes))

