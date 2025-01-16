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

freq_impute <- function(gl, mt, ploidity){
  mt <- as.matrix(gl)
  # Get the allelic frequencies
  q_allele <- adegenet::glMean(gl)
  # Linear index of nas
  idx_na <- which(is.na(mt))
  na_loc_idx <- sapply(idx_na, function(x){
    loc_idx <- ceiling(x/nrow(mt))
    return(loc_idx)
  })
  
  imp <- unname(unlist(lapply(q_allele[na_loc_idx],
                              function(x) {return(as.numeric(i_freq_impute(q_frq = x, ploidity)))})))
  return(split(idx_na, imp))
}

apply_imputation <- function(mt, imp_dict){
  for (dosage in names(imp_dict)) {
    # Convert the list name to a numeric value
    num_dosage <- as.numeric(dosage)
    
    # Get the linear indices associated with this value
    idx <- imp_dict[[dosage]]
    
    # Assign the value to these positions in the matrix
    mt[idx] <- num_dosage
  }
  return(mt)
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
  
  loci_all_nas <- sum(adegenet::glNA(gl)/ploidity == adegenet::nInd(gl))
  nas_number <- sum(adegenet::glNA(gl))/ploidity
  number_imputations <- nas_number - (loci_all_nas * adegenet::nInd(gl))
  
  mt <- as.matrix(gl)
  
  if(loci_all_nas > 0){
    cli::cli_warn("There are {loci_all_nas} loci with all missing data")
  }
  
  cli::cli_inform("Missing genotype calls {number_imputations}")
  
  if(method == 'frequency'){
    imp_dict <- freq_impute(gl, mt, ploidity)
  }
  
  # apply the imputation creating a new gl instance
  imp_mt <- apply_imputation(mt, imp_dict)
  
  
  imp_gl <- new("genlight",
            imp_mt,
            ploidy = ploidity,
            loc.names = gl@loc.names,
            ind.names = gl@ind.names,
            chromosome = gl@chromosome,
            position = gl@position)
  
  adegenet::alleles(imp_gl) <- adegenet::alleles(gl)
  imp_gl <- recalc_metrics(imp_gl)
  
  return(list(gl = imp_gl, log = imp_dict))
}
