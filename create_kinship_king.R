# system("~/bin/king/king -b ./data/plink_files/atchely_imputed.bed --kinship --prefix atchely_imputed")
#
# within_fam = read_tsv("atchely_imputed.kin", col_types = "dddidddddd")
# between_fam = read_tsv("atchely_imputed.kin0", col_types = "ddddiddd")
#
# within_fam$ID1 = as.character(within_fam$ID1)
# within_fam$ID2 = as.character(within_fam$ID2)
#
# between_fam$ID1 = as.character(between_fam$ID1)
# between_fam$ID2 = as.character(between_fam$ID2)
#
# snpKin = kinship(pedAll)
# for(i in 1:nrow(within_fam)){
#   snpKin[within_fam$ID1[i], within_fam$ID2[i]] = within_fam$Kinship[i]
# }
# for(i in 1:nrow(between_fam)){
#   snpKin[between_fam$ID1[i], between_fam$ID2[i]] = between_fam$Kinship[i]
# }
#
# write_tsv(tbl_df(snpKin), "./data/king_snp_kinship_matrix.tsv")
snpKin = as.matrix(read_tsv("./data/king_snp_kinship_matrix.tsv"))
dimnames(snpKin)[1] = dimnames(snpKin)[2]
