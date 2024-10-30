
# filter by missing rate by sample a geneind obj
filter_missing_rate_by_indv <- function(gl, threshold){
  index <- gl@other$ind.metrics$ind_miss <= threshold
  gl2 <- gl[index,]
  gl2 <- recalc_metrics(gl2)
  return(gl2)
}

# filter by missing rate by marker a geneind obj
filter_missing_rate_by_marker <- function(gl, threshold){
  index <- gl@other$loc.metrics$loc_miss <= threshold
  gl2 <- gl[, index]
  gl2 <- recalc_metrics(gl2)
  return(gl2)
}

filter_MAF <- function(gl, threshold){
  index <- gl@other$loc.metrics$maf >= threshold
  gl2 <- gl[, index]
  gl2 <- recalc_metrics(gl2)
}

filter_heterozygosis_by_loc <- function(gl, threshold){
  index <- gl@other$loc.metrics$loc_het <= threshold
  gl2 <- gl[, index]
  gl2 <- recalc_metrics(gl2)
}

filter_heterozygosis_by_ind <- function(gl, threshold){
  index <- gl@other$ind.metrics$ind_het <= threshold
  gl2 <- gl[index,]
  gl2 <- recalc_metrics(gl2)
}



