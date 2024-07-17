
# Missingness metrics -----------------------------------------------------

#' Get Locus Missingness
#'
#' Compute the missing rate for each locus (marker) in a genlight object.
#' The missing rate is a value between 0 and 1, where 0 indicates no missing data
#' for that locus, and 1 indicates that all samples have missing data for that locus.
#'
#' @param gl A genlight object.
#'
#' @return A numeric vector of length equal to the number of loci (markers),
#'   containing the missing rate for each locus.
#'
#' @export
#'
#' @examples
#' data(example_genlight)
#' loc_miss <- get_loc_missing(example_genlight)
#' head(loc_miss)
get_loc_missing <- function(gl) {
  # Get the number of occurrences of NAs (missing data) for each marker
  NA_counts <- adegenet::glNA(gl)
  
  # Divide the NA counts by the total number of samples to get the missing rate
  NA_counts <- NA_counts / adegenet::nInd(gl)/max(adegenet::ploidy(gl))
  
  return(NA_counts)
}

#' Get Individual Missingness
#'
#' Compute the missing rate for each individual (sample) in a genlight object.
#' The missing rate is a value between 0 and 1, where 0 indicates no missing data
#' for that individual, and 1 indicates that all loci have missing data for that individual.
#'
#' @param gl A genlight object.
#'
#' @return A numeric vector of length equal to the number of individuals (samples),
#'   containing the missing rate for each individual.
#'
#' @export
#'
#' @examples
#' data(example_genlight)
#' ind_miss <- get_ind_missing(example_genlight)
#' head(ind_miss)
get_ind_missing <- function(gl) {
  # Convert the genlight object to a matrix and identify missing genotype calls
  mt <- is.na(as.matrix(gl))
  
  # Calculate the proportion of missing data for each individual (row)
  ind_miss <- rowSums(mt) / adegenet::nLoc(gl)
  
  return(ind_miss)
}


#' Get Overall Missingness
#'
#' Compute the overall missing rate for a genlight object.
#' The overall missing rate is the proportion of missing genotype calls
#' across all individuals and loci in the dataset.
#'
#' @param gl A genlight object.
#'
#' @return A single numeric value representing the overall missing rate.
#'
#' @export
#'
#' @examples
#' data(example_genlight)
#' overall_miss <- get_overall_missingness(example_genlight)
#' print(overall_miss)
get_overall_missingness <- function(gl) {
  # Convert the genlight object to a matrix
  mt <- as.matrix(gl)
  
  # Identify missing genotype calls
  mt <- is.na(mt)
  
  # Calculate the overall missing rate
  overall_miss <- sum(mt) / (nrow(mt) * ncol(mt))
  
  return(overall_miss)
}


# Heterozygosity metrics --------------------------------------------------

get_heterozygosity_metrics <- function(gl, ploidy = 2){
  # Boolean matrix of genotype calls where 0 > dosage < ploidy
  mt <- as.matrix(gl)
  het_ind_loc <- mt > 0 & mt < ploidy
  het_loc <- colSums(het_ind_loc , na.rm = T)/adegenet::nInd(gl)
  het_ind <- rowSums(het_ind_loc , na.rm = T)/adegenet::nLoc(gl)
  return(list(het_ind = het_ind, het_loc = het_loc))
}

get_loc_heterozygosity <- function(gl, ploidy = 2){
  mt <- as.matrix(gl)
  het_ind_loc <- mt > 0 & mt < ploidy
  het_loc <- colSums(het_ind_loc , na.rm = T)/adegenet::nInd(gl)
  return(het_loc)
}

get_ind_heterozygosity <- function(gl, ploidy = 2){
  mt <- as.matrix(gl)
  het_ind_loc <- mt > 0 & mt < ploidy
  het_ind <- rowSums(het_ind_loc , na.rm = T)/adegenet::nLoc(gl)
  return(het_ind)
}


get_maf <- function(gl){
  maf <- adegenet::glMean(gl)*(1/2)
  return(maf)
}

recalc_metrics <- function(gl){

  gl@other$loc.metrics <- data.frame(
    maf = get_maf(gl),
    loc_miss = get_loc_missing(gl),
    loc_het = get_loc_heterozygosity(gl, max(adegenet::ploidy(gl)))
  )
  
  
  gl@other$ind.metrics <- data.frame(
                                     ind_miss = get_ind_missing(gl),
                                     ind_het = get_ind_heterozygosity(gl, max(adegenet::ploidy(gl)))
                                     )
  return(gl)
}




