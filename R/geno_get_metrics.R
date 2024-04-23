get_freq_het <- function(geno_obj, plodity = 2){
  # Get allelic dosage
  d <- allelic_dosage(geno_obj$geno)
  freq_het <- colMeans(x > 0 & d < plodity , na.rm = T)
  return(freq_het)
}

get_loc_missing <- function(geno_obj){
  # Get allelic dosage
  d <- allelic_dosage(geno_obj$geno)
  na_count <- colSums(is.na(d), na.rm = T)
  return(na_count/adegenet::nInd(geno_obj$geno))
}

get_maf <- function(geno_obj){
  # Get allelic dosage
  d <- allelic_dosage(geno_obj$geno)
  alf <- dartR::gl.alf(d)[, 2]
  maf <- ifelse(alf > 0.5, 1 - alf, alf)
  return(maf)
}
