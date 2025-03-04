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


#' Mask a fraction of genotype calls for measure the 
#' accuracy of an imputation method
#'
#' @param gl A genlight instance
#' @param fraction percentage of data to mask [0 - 1], default = 0.1
#'
#' @return A genlight instance where the genotype matrix have as NA
#' the given fraction of missing data
#' @export
#'
#' @examples
set_fraction_na <- function(gl, fraction = 0.1) {
    # Ensure fraction is between 0 and 1
    if (!is.numeric(fraction) || fraction < 0 || fraction > 1) {
      stop("fraction must be a number between 0 and 1")
    }
    m <- as.matrix(gl)
    
    # Get total number of elements
    total_elems <- length(m)
    
    # Compute how many elements to replace
    n_replace <- round(total_elems * fraction)
    
    # Randomly choose positions
    positions <- sample(seq_len(total_elems), size = n_replace, replace = FALSE)
    
    # Replace chosen positions with NA
    m[positions] <- NA
    
    masked_gl <- new("genlight", m, ploidy=max(gl@ploidy))
    masked_gl <- recalc_metrics(masked_gl)
    return(masked_gl)
}


#' Compare the genotype calls (gt) of a reference genlight object
#' against an imputed gl. Measure the accuracy of the imputation
#' process comparing the reference gt with imputed gt
#'
#' @param ref_gl Rreference genlight instance
#' @param imp_gl Imputed genlight instance
#'
#' @return
#' @export
#'
#' @examples
get_accuracy <- function(ref_gl, imp_gl){
  # Get the gt matix of the two gls
  ref_m <- as.matrix(ref_gl)
  imp_m <- as.matrix(imp_gl)
  # Boolean matrix where the elements match/missmatch
  comp_m <- ref_m == imp_m
  # Boolean matrix where the element is NA
  na_m <- is.na(ref_m)
  target_cells <- which(na_m)
  
  # Convert to NA the element where in ref is NA (not comparable)
  comp_m[target_cells] <- NA
  # Sum by column the correctly imputed gt
  matchs_sums <- colSums(comp_m, na.rm = T)
  accuracy <- matchs_sums/colSums(!na_m)
  
  return(accuracy)
}



