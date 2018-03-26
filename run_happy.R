detach("package:happy", unload = TRUE)
library(happy)

h <- happy( './data/f5f6_genotypes_happy_chr6.PED', 'data/markers_strain_split_chr6.txt',
            generations=6,
            file.format = "ped",
            phase ="estimate" )
qtl_fit = hfit(h, model = "additive")
happyplot(qtl_fit)

str(h$map)
h$markers[1:length(h$markers)-1]
x = hdesign( h, h$markers[1], model='full' )
zapsmall(x)
x[1500,]
ind1 = laply(1:((length(h$markers))-1),
      function(i) hdesign( h, i, model='additive' )[900,])
gather(data.frame(ID = 1:dim(ind1)[1], ind1), key, value, -ID) %>%
  ggplot(aes(ID, value, color = key)) + geom_line()
zapsmall(hdesign( h, 1, model='additive' ), 3)
happy:::strain.effects(h, qtl_fit)

h$subjects
happy.matrices(h)

h$nam
h$haploid
is.null(h$matrices)

zapsmall(hdesign(h, 3), 3)
zapsmall(hprob(h, 1))
