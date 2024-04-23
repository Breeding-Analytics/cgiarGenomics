read_tabular_geno <- function(path, sep = '\t'){
  # missing data representation values
  missingData = c("NN","FAIL","FAILED","Uncallable","Unused","-","NA","",-9)

  df <- as.data.frame(data.table::fread(path,
                                        sep = sep,
                                        header = TRUE,
                                        na.strings = missingData))
  return(df)
}

print_log_message <- function(message, log_level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  formatted_message <- paste("[", timestamp, "][", log_level, "]: ", message, sep = "")
  cat(formatted_message, "\n")
}

addUnits <-
  function(n) {
    labels <-
      ifelse(n < 1000, n,  # less than thousands
             ifelse(n < 1e6, paste0(round(n/1e3), 'k'),  # in thousands
                    ifelse(n < 1e9, paste0(round(n/1e6), 'M'),  # in millions
                           ifelse(n < 1e12, paste0(round(n/1e9), 'B'), # in billions
                                  'too big!'))))
    return(labels)
  }

get_midpoint <- function(cut_label) {
  mean(as.numeric(unlist(strsplit(gsub("\\(|\\)|\\[|\\]", "", as.character(cut_label)), ","))))
}

allelic_dosage <- function(genind){
  # convert to genind data structure using a separator empty
  locna <- genind@loc.n.all
  ccc <- 1

  # keep the less frequent allele countings
  for (i in 2:length(locna)) {
    if (locna[i - 1] == 1) {
      ccc[i] <- ccc[i - 1] + 1
    }
    else {
      ccc[i] <- ccc[i - 1] + 2
    }
  }

  alellelic_dosage <- genind@tab[, ccc]
  return(alellelic_dosage)
}








get_file_size <- function(path){
  file_info <- file.info(path)
  # Extract file size in bytes
  file_size_bytes <- file_info$size
  return(file_size_bytes)
}



get_overall_missingness <- function(gdsobj, samples_id, snps_id){
  mt <- snpgdsGetGeno(gdsob,
                      sample.id = samples_id,
                      snp.id = snps_id)
  total_gt <- dim(mt)[1]*dim(mt)[2]
  missing_calls <- sum(is.na(SNPRelate::snpgdsGetGeno(genofile)))
  return(missing_calls/total_gt)
}



sample_gt_mt <- function(gt, fraction = 0.1){
  non_na_indices <- which(!is.na(gt))
  total_non_na_elements <- length(non_na_indices)
  elements_to_sample <- round(fraction * total_non_na_elements)
  if(elements_to_sample > 10e3 ){
    elements_to_sample <- 10e3
  }
  # Randomly sample the specified percentage of non-NA elements
  sampled_indices <- sample(non_na_indices, size = elements_to_sample)
  return(gt[sampled_indices])
}

validate_gt_numeric <- function(lst) {
  result <- NULL
  for (element in lst) {
    if (grepl("^-?\\d+$", element)) {
      result <- c(result, TRUE)
    } else {
      result <- c(result, FALSE)
    }
  }
  return(all(result))
}

validate_ploidity <- function(lst, ploidity = 2){
  if(validate_gt_numeric(lst)){
    # is numeric
    if(max(lst) != ploidity){
      return(FALSE)} else {
        return(TRUE)
      }
  } else {
    # is character
    lengths_nn <- unlist(lapply(lst, nchar))
    if(all(lengths_nn == ploidity)){
      return(TRUE)} else {
        return(FALSE)
      }
  }
}


validate_hapmap <- function(table){

  hapmap_snp_attr <- c('rs#', 'alleles', 'chrom', 'pos', 'strand', 'assembly#',
                       'center', 'protLSID', 'assayLSID', 'panelLSID', 'QCcode',
                       'rs', 'assembly','panel')

  validator <- list(format = "hapmap")

  if(length(intersect(hapmap_snp_attr, colnames(table)[1:11])) != 11){
    validator$validated <- FALSE
    validator$error <- 'Column names doesn\'t match to hapmap format'
  }

  return(validator)
}
