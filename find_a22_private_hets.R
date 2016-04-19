source("find_het_sites.R")

A22_het_pos = just_snps %>% filter(A22_hets == 1)

A22_hets_array_list = dlply(A22_het_pos, .(CHROM), function(x) as.matrix(select(x, A13_1:A42_2)))

x = A22_hets_array_list[[1]][10,]
A22_alele1_mask = as.character((names(x) == "A22_1") + 1)
A22_alele2_mask = (names(x) == "A22_2") + 1

get_private_het = function(x){
       if(all(tapply(x, A22_alele1_mask, function(row) length(unique(row))) == c(1, 1)))
              return(1)
  else if(all(tapply(x, A22_alele2_mask, function(row) length(unique(row))) == c(1, 1)))
          return(1)
    else
            return(0)
}
apply(A22_hets_array_list[[1]], 1, get_private_het)

is_private_heterozygote <-
  ldply(A22_hets_array_list,
        function(A22_het_array) adply(A22_het_array,
                                      1, get_private_het),
        .parallel = TRUE)
