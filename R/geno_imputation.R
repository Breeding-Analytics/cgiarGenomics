#' Imputation with allele frequency
#' 
#' Assuming a bi-allelic marker, using the observed allelic frequency for one allele
#' is sampled the genotype call for any ploidity level
#'
#' @param q_frq 
#' @param ploidity 
#'
#' @return
#' @export
#'
#' @examples
i_freq_impute <- function(q_frq, ploidity = 2){
  if(!rlang::is_integerish(ploidity)){
    cli::cli_abort("`ploidity` is not an integer: {ploidity}")  
  }
  dosage <- 0
  for(i_chromatid in seq(ploidity)){
    i_dosage <- sample(c(1,0), size = 1, 
                       prob = c(q_frq, 1 - q_frq), replace = T)
    dosage <- dosage + i_dosage
  }
  return(dosage)
}


#' Impute a gl object
#'  
#' Impute a genlight object using currently the frequency method.
#' Is returned a modification dataframe with the sample, locus id and the imputed
#' genotype call
#'
#' @param gl 
#' @param ploidity 
#' @param method 
#'
#' @return
#' @export
#'
#' @examples
impute_gl <- function(gl, ploidity = 2, method = 'frequency'){
  
  loci_all_nas <- sum(adegenet::glNA(gl)/ploidity > adegenet::nInd(gl))
  nas_number <- sum(adegenet::glNA(gl))/ploidity
  number_imputations <- nas_number - (loci_all_nas * adegenet::nInd(gl))
  
  if(loci_all_nas > 0){
    cli::cli_warn("There are {loci_all_nas} loci with all missing data")
  }
  
  cli::cli_inform("Missing genotype calls {number_imputations}")
  
  mt <- as.matrix(gl)
  q_allele <- adegenet::glMean(gl)
  
  loc_na <- as.data.frame(which(is.na(mt), arr.ind = TRUE, useNames = F))
  
  
  if(method == 'frequency'){
    imp <- unname(unlist(lapply(q_allele[loc_na[,2]],
                                function(x) {return(as.numeric(i_freq_impute(q_frq = x, ploidity)))})))
  }
  colnames(loc_na) <- c('row', 'col')
  loc_na$call <- imp
  return(loc_na)
}
