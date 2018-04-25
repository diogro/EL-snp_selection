if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}

positions = read.table("~/projects/EL-snp_selection/data/marker_positions.txt", header = FALSE, stringsAsFactors = FALSE) %>% tbl_df
names(positions) = c("chr2", "pos", "chr", "gpos")

plink_map_raw = read_table2("~/projects/EL-snp_selection/data/plink_files/imputed_autossomal.bim",
                            comment = "#", col_names = c("chr", "ID", "gpos", "pos", "minor", "major"))

plink_map =
  inner_join(select(plink_map_raw, -gpos),
             positions,
             by = c("chr", "pos")) %>%
  mutate(gpos = gpos) %>%
  dplyr::select(chr, ID, gpos, pos, minor, major)

all(plink_map[,c(1, 2, 4, 5, 6)] == plink_map_raw[,c(1, 2, 4, 5, 6)])

write_delim(plink_map, "~/projects/EL-snp_selection/data/plink_files/imputed_autossomal.bim", col_names = FALSE, delim = "\t")