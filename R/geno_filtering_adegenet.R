# filter monomorphic markers
filter_monomorphic <- function(geno_obj)  {

  gi <- geno_obj$geno
  # select poly markers
  polymorphic <- names(which(adegenet::isPoly(gi)))
  # number of total loci
  nloci <- adegenet::nLoc(gi)
  if (length(polymorphic) > 0) {
    # creates a new geneind obj with only poly markers
    poly_gi <- gi[loc=polymorphic]
    msj <- paste("From", nloci,
          "loci were removed",
          nloci - length(polymorphic),
          "monomorphic Loci!", sep = ' ')
    print_log_message(msj)

    id_var <- "ID"
    meta <- geno_obj$meta %>%
      filter(!!id_var %in% names(poly_gi@all.names))

    geno_obj <- list(geno = poly_gi, meta = meta)

    return(geno_obj)
  } else{
    print_log_message("Any monomorphic marker was found!")
    return(geno_obj)
  }
}

# filter by missing rate by marker a geneind obj
filter_missing_rate_by_marker <- function(geno_obj, threshold){
  gi <- geno_obj$geno
  lmiss_f_gi <- poppr::missingno(gi,
                   type = 'loci',
                   cutoff = threshold)
  id_var <- "ID"
  meta <- geno_obj$meta %>%
    filter(!!id_var %in% names(lmiss_f_gi@all.names))

  geno_obj <- list(geno = lmiss_f_gi, meta = meta)

  return(geno_obj)
}


# filter by missing rate by sample a geneind obj
filter_missing_rate_by_indv <- function(geno_obj, threshold){
  gi <- geno_obj$geno
  smiss_f_gi <- poppr::missingno(gi,
                                 type = 'geno',
                                 cutoff = threshold)
  id_var <- "ID"
  meta <- geno_obj$meta %>%
    filter(!!id_var %in% names(smiss_f_gi@all.names))

  geno_obj <- list(geno = smiss_f_gi, meta = meta)

  return(geno_obj)
}


# filter by minor allele frequency
filter_MAF <- function(geno_obj, threshold){
  gi <- geno_obj$geno
  maf_loci <- names(which(adegenet::minorAllele(gi) >= threshold))
  # number of total loci
  nloci <- adegenet::nLoc(gi)
  if (length(maf_loci) > 0) {
    # creates a new geneind obj with maf greater than threshold loci
    maf_f_gi <- gi[loc = maf_loci]
    msj <- paste("From", nloci,
                 "loci were removed",
                 nloci - length(maf_loci),
                 "Loci because MAF <",
                 threshold,"!", sep = ' ')
    print_log_message(msj)
    id_var <- "ID"
    meta <- geno_obj$meta %>%
      filter(!!id_var %in% names(maf_f_gi@all.names))

    geno_obj <- list(geno = maf_f_gi, meta = meta)

    return(geno_obj)}
  else{
    msg <- paste("Any locus with MAF smaller than",
                 threshold,
                 " was found!", sep = ' ')
    print_log_message(msj)
    return(geno_obj)
  }
}


# Adegenet native implementation for filter
# # filter by missing rate by marker a geneind obj
# filter_missing_rate_by_marker <- function(gi, threshold){
#   miss_samples <- names(which(adegenet::propTyped(gi, by = 'ind') <= threshold))
#   # number of total samples
#   nSamp <- adegenet::nInd(gi)
#   if (length(miss_samples) > 0) {
#     # creates a new geneind obj with samples with miss rate smaller than threshold
#     smiss_f_gi <- gi[miss_samples,]
#     msg <- paste("From", nSamp,
#                  "samples were removed",
#                  nSamp - length(smiss_f_gi),
#                  "samples because missing rate >",
#                  threshold,"!", sep = ' ')
#     print_log_message(msj)
#     return(smiss_f_gi)}
#   else{
#     msg <- paste("Any sample with missing rate >",
#                  threshold," was found!", sep = ' ')
#     print_log_message(msg)
#     return(gi)
#   }
# }
# system.time(filter_missing_rate_by_indv(gi, 0.1))
#


