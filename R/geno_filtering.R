filter_missing_rate_by_indv <- function(gdsobj, threshold, samples_id = NULL, snps_id = NULL){
  missing_rate <- SNPRelate::snpgdsSampMissRate(gdsobj,
                                snp.id = snps_id,
                                sample.id = samples_id,
                                with.id = TRUE)
  filter <- missing_rate <= threshold
  retained_set <- names(filter)[filter]

  return(retained_set)
}

filter_missing_rate_by_marker <- function(gdsobj, threshold, samples_id = NULL, snps_id = NULL){
  missing_rate <- SNPRelate::snpgdsSNPRateFreq(gdsobj,
                                                snp.id = snps_id,
                                                sample.id = samples_id,
                                                with.snp.id = TRUE,
                                                with.id = TRUE)
  filter <- missing_rate$MissingRate <= threshold
  retained_set <- missing_rate$snp.id[filter]
  return(retained_set)
}
