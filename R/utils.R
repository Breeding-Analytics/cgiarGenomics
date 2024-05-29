#' Read a text table and use common na representations for convert to NA
#'
#' @param path Path to tabular data
#' @param sep  Delimitator of the table cells
#'
#' @return Returns a dataframe
#' @export
#'
#' @examples
read_tabular_geno <- function(path, sep = '\t'){
  # missing data representation values
  missingData = c("N","NN","FAIL","FAILED","Uncallable","Unused","-","NA","",-9)

  df <- as.data.frame(data.table::fread(path,
                                        sep = sep,
                                        header = TRUE,
                                        na.strings = missingData))
  return(df)
}

#' Get allelic counts for a single variant according to ploidity level
#' INDELS are not allowed
#'
#' @param v Vector containing the genotype calls for each sample in character format
#' @param ploidity Ploitidy level
#' @param sep Separator used in genotype call for each chromatid
#'
#' @return Return a table with allelic counts
#' @export
#'
#' @examples
#' v <- c("A/A/A/A", "A/C/A/A","C/C/A/A","C/C/C/C")
#' get_alleles_count_char(v, ploidity = 4, sep = "/")
get_alleles_count_char <- function(v, ploidity = 2,  sep = ""){
  vu <- v
  # Removes all missing genotype calls
  v <- v[which(!is.na(v))]

  # Remove the separator from the genotype call A/A -> AA
  v <- paste(v, collapse = " ")
  v <- gsub(sep, "", v)
  v <- unlist(strsplit(v, " "))

  # Count by genotype call class
  tmp <- table(v)
  tmp <- tmp[which(!is.na(tmp))]
  geno_class <- names(tmp)

  # verify ploidity if the length of genotype classes are divisible by
  # Ploidity AA % 2 = 0  or AABB % 4 = 0
  poly_check <- sapply(names(tmp), nchar) %% ploidity


  if(sum(poly_check > 0 )){
    print(tmp)
    print(vu)
    stop(paste0("Marker with polodity different to: ", ploidity))
  } else {
    # Get the unique set of alleles dividing the genotype classes according
    # to ploidity AABB and ploidity = 4 => A,A,B,B
    alleles_mt <- split_strings(names(tmp), ploidity)
    alleles <- unique(as.vector(alleles_mt))
    # Initialize he counts list
    counts <- rep(0, length(alleles))

    # Iterate in each genotype class and update the alleles counts in
    # the counts list
    for(gclass in 1:ncol(alleles_mt)){
      ialleles <- alleles_mt[,gclass]
      for(iallele_idx in 1:length(ialleles)){
        counts[which(alleles == ialleles[iallele_idx])] <- counts[which(alleles == ialleles[iallele_idx])] + as.numeric(tmp[gclass])
      }
    }
    # Set names to alleles vector
    names(counts) <- alleles
    # Sort allele counts
    counts <- counts[order(as.numeric(counts), decreasing = T)]
    return(counts)
  }
}

print_log_message <- function(message, log_level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  formatted_message <- paste("[", timestamp, "][", log_level, "]: ", message, sep = "")
  cat(formatted_message, "\n")
}

#' Split a genotype call strings according to the ploidity level ex:
#' Ploidity call alleles. Indels are not allowed
#' 2  c(AA, BB) => A, B
#'                 A, B
#'
#' 4 c(AAAA,AABB) =>  A, A
#'                    A, A
#'                    A, B
#'                    A, B
#'
#' @param strings Vector of genotype calls
#' @param ploidity Level of ploidity
#'
#' @return A matrix with columns each genotype call and rows alleles
#' @export
#'
#' @examples
split_strings <- function(strings, ploidity = 2 ) {
  sapply(strings, function(s, ploidity) {
    if(!is.na(nchar(s))){
      split_length <- nchar(s)/ploidity
      substring(s, seq(1, nchar(s), split_length), seq(split_length, nchar(s), split_length))
    }}, ploidity = ploidity)
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
