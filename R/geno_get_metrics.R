get_loc_het <- function(geno_obj, plodity = 2){
  # select the lowest value, in this case should be the same for both alleles
  f1 <- function(vec) {
    vec <- sort(vec, decreasing = TRUE)
    res <- vec[min(2, length(vec))]
    return(res)
  }

  het_ind_loc <- adegenet::tab(geno_obj$geno) > 0 & adegenet::tab(geno_obj$geno) < plodity
  freq_het <- colMeans(het_ind_loc , na.rm = T)
  out <- tapply(freq_het, adegenet::locFac(geno_obj$geno), f1)
  return(out)
}

get_ind_het <- function(geno_obj, plodity = 2){
  het_ind_loc <- adegenet::tab(geno_obj$geno) > 0 & adegenet::tab(geno_obj$geno) < plodity
  out <- rowMeans(het_ind_loc , na.rm = T)
  return(out)
}

get_loc_missing <- function(geno_obj){
  # select the lowest value, in this case should be the same for both alleles
  f1 <- function(vec) {
    vec <- sort(vec, decreasing = TRUE)
    res <- vec[min(2, length(vec))]
    return(res)
  }
  na_count <- colSums(is.na(adegenet::tab(geno_obj$geno)), na.rm = T)
  out <- tapply(na_count, adegenet::locFac(geno_obj$geno), f1)/adegenet::nInd(geno_obj$geno)
  return(out)
}

get_ind_missing <- function(geno_obj){
  # select the lowest value, in this case should be the same for both alleles
  f1 <- function(vec) {
    vec <- sort(vec, decreasing = TRUE)
    res <- vec[min(2, length(vec))]
    return(res)
  }
  na_count <- rowSums(is.na(adegenet::tab(geno_obj$geno)), na.rm = T)
  out <- na_count/adegenet::nLoc(geno_obj$geno)
  return(out)
}

get_maf <- function(geno_obj){

  f1 <- function(vec) {
    vec <- sort(vec, decreasing = TRUE)
    res <- vec[min(2, length(vec))]
    return(res)
  }

  freq <- apply(adegenet::tab(geno_obj$geno, freq = TRUE), 2, mean, na.rm = TRUE)
  out <- tapply(freq, adegenet::locFac(geno_obj$geno), f1)
  return(out)
}



