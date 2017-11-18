#===== Prepare some data

#===== Test calculate group size

calc_group_size(N = 270, donor2cancer.use$histology_abbreviation, is.balance = T)
calc_group_size(N = 271, donor2cancer.use$histology_abbreviation, is.balance = T)
calc_group_size(N = 100, donor2cancer.use$histology_abbreviation, is.balance = T)
calc_group_size(N = 270, donor2cancer.use$histology_abbreviation, is.balance = F)
calc_group_size(N = 1188, donor2cancer.use$histology_abbreviation, is.balance = F)


#===== Test sampling hits matrix
library(Matrix)
library(profvis)
ase <- as.matrix(ase)
profvis(
aa <- sampling_hits_matrix(hits.mat = ase, n = 100, n.repeat = 10, group.table = donor2cancer.use, n.cores = 2, do.parallel = F)
)
pryr::object_size(aa)
sampling_by_group
library(do)

remove.packages("doMPI")