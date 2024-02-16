get_overall_missingness <- function(gdsobj, samples_id, snps_id){
  mt <- snpgdsGetGeno(gdsob,
                sample.id = samples_id,
                snp.id = snps_id)
  total_gt <- dim(mt)[1]*dim(mt)[2]
  missing_calls <- sum(is.na(SNPRelate::snpgdsGetGeno(genofile)))
  return(missing_calls/total_gt)
}
