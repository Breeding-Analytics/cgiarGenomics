
# Purity Metrics ----------------------------------------------------------
get_purity <- function(gl, inds) {
  if (sum(inds %in% indNames(gl)) < length(inds)) {
   cli::cli_abort("At least one duplicated ind not in gl")
  }
  ind_idx <- which(inds %in% indNames(gl))
  cmb <- combn(inds, 2)
  ind_x  <- cmb[1, ]
  ind_y <- cmb[2, ]
  comparision <- purrr::map2_dfr(ind_x, ind_y, ~{
    ind_idx <- which(indNames(gl) %in% c(.x,.y))
    pair_gl <- gl[ind_idx,]
    paired_ind_comparision(pair_gl,.x,.y)
  })
  return(comparision)
}

ibs_matrix_purrr <- function(gl, ploidy) {
  G <- t(as.matrix(gl))
  cols <- as.data.frame(G, check.names = FALSE)
  ids  <- colnames(cols)
  
  score_list <- map(cols, function(x)
    map_dbl(cols, function(y) ibs_dosage(x, y, ploidy)$score)
  )
  overlap_list <- map(cols, function(x)
    map_int(cols, function(y) ibs_dosage(x, y, ploidy)$overlap)
  )
  
  score_mat   <- do.call(cbind, score_list)
  overlap_mat <- do.call(cbind, overlap_list)
  
  rownames(score_mat) <- ids; colnames(score_mat) <- ids
  rownames(overlap_mat) <- ids; colnames(overlap_mat) <- ids
  
  # numeric safety for diagonal (optional)
  diag(score_mat) <- 1
  
  return(list(score = score_mat, overlap = overlap_mat))
}

ibs_dosage <- function(x, y, ploidy) {
  if (length(x) != length(y)) stop("x and y must have the same length")
  if (length(ploidy) != 1 || !is.numeric(ploidy) || ploidy <= 0) {
    stop("ploidy must be a positive number")
  }
  # mask non-missing pairs
  ok <- !is.na(x) & !is.na(y)
  n_ok <- sum(ok)
  if (n_ok == 0L) return(list(score = NA, overlap = 0L))
  
  xv <- x[ok]; yv <- y[ok]
  
  # locus-wise IBS: 1 - |Δdosage| / P
  ibs_locus <- 1 - abs(xv - yv) / ploidy
  # numeric safety clamp
  ibs_locus <- pmin(pmax(ibs_locus, 0), 1)

  list(score = mean(ibs_locus), overlap = n_ok)
}

get_paired_IBS <- function(gl, ploidy) {
  ibs <- ibs_matrix_purrr(gl, ploidy)
  return(ibs)
}

paired_ind_comparision <- function(gl, x, y) {
  if (length(indNames(gl)) > 2) {
    cli::cli_abort("gl have more than two inds, not paired")
  }
  inf_loc_idx <- unname(which(get_loc_missing(gl) == 0))
  tgt_gl <- gl[,inf_loc_idx]
  mt <- as.matrix(tgt_gl)
  ploidity_lvl <- max(adegenet::ploidy(gl))
  hom_comp_idx <- which(colSums(mt) %in% c(0, 2*ploidity_lvl))
  het_comp_idx <- which(!colSums(mt) %in% c(0, 2*ploidity_lvl))
  
  hom_diff <- sum(mt[1,hom_comp_idx] != mt[2,hom_comp_idx])
  het_diff <- sum(mt[1,het_comp_idx] != mt[2,het_comp_idx])
  out <- data.frame(ind1 = x,
                    ind2 = y,
              sites = dim(mt)[2],
              hom_diff = hom_diff,
              het_diff = het_diff)
  return(out)
}


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
  ind_miss <- Matrix::rowSums(mt) / adegenet::nLoc(gl)
  
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
  het_loc <- Matrix::colSums(het_ind_loc , na.rm = T)/adegenet::nInd(gl)
  het_ind <- Matrix::rowSums(het_ind_loc , na.rm = T)/adegenet::nLoc(gl)
  return(list(het_ind = het_ind, het_loc = het_loc))
}

get_loc_heterozygosity <- function(gl, ploidy = 2){
  mt <- as.matrix(gl)
  het_ind_loc <- mt > 0 & mt < ploidy
  het_loc <- Matrix::colSums(het_ind_loc , na.rm = T)/adegenet::nInd(gl)
  return(het_loc)
}

get_ind_heterozygosity <- function(gl, ploidy = 2){
  mt <- as.matrix(gl)
  het_ind_loc <- mt > 0 & mt < ploidy
  het_ind <- Matrix::rowSums(het_ind_loc , na.rm = T)/adegenet::nLoc(gl)
  return(het_ind)
}


get_maf <- function(gl, ploidy = 2){
  alf <- adegenet::glMean(gl)*(1/ploidy)
  maf <- ifelse(alf > 0.5, 1 - alf, alf)
  return(maf)
}

get_inbreeding <- function(gl){
  p <- 1 - gl@other$loc.metrics$maf
  he <- 2 * gl@other$loc.metrics$maf * p
  Fis <- 1 - (gl@other$loc.metrics$loc_het/he)
  return(Fis)
}

recalc_metrics <- function(gl){

  gl@other$loc.metrics <- data.frame(
    maf = get_maf(gl, max(adegenet::ploidy(gl))),
    loc_miss = get_loc_missing(gl),
    loc_het = get_loc_heterozygosity(gl, max(adegenet::ploidy(gl)))
  )
  
  # Use already calculated loc stats
  gl@other$loc.metrics$loc_Fis <- get_inbreeding(gl)
  
  gl@other$ind.metrics <- data.frame(
                                     ind_miss = get_ind_missing(gl),
                                     ind_het = get_ind_heterozygosity(gl, max(adegenet::ploidy(gl)))
                                     )
  return(gl)
}


get_overall_summary <- function(gl){
  ninds <- adegenet::nInd(gl)
  nlocs <- adegenet::nLoc(gl)
  overall_missiness <- mean(gl@other$ind.metrics$ind_miss)
  overall_heterozygosity <- mean(gl@other$ind.metrics$ind_het)
  overall_maf <- mean(gl@other$loc.metrics$maf, na.rm = T)
  
  out <- list(
    nind = ninds,
    nloc = nlocs,
    ov_miss = overall_missiness,
    ov_het = overall_heterozygosity,
    ov_maf = overall_maf
  )
  return(out)
}



