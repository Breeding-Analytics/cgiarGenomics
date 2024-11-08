#' @param geno_metadata Dataframe with the expected columns id, chrom, pos, ref, alt
#'
#' @return Same input dataframe but added a `filter` column if true the marker don't meet
#' the minimal requirements: no duplicated markers, physical position, bi-allelic, 
#'  allelic information, SNP.
#' @export
#'
#' @examples
#' geno_metadata(meta_daya)
process_metadata <- function(geno_metadata){
  # Expected fields id Chrom pos ref alt
  meta_columns <- c('id', 'chrom', 'pos', 'ref', 'alt')
  
  # Check if the columns in the input file match the expected column names
  if (length(intersect(meta_columns, colnames(geno_metadata))) != 5) {
    cli::cli_abort("metadata file doesn't have the expected column names,
                   have: {colnames(geno_metadata)}")
  }
  
  geno_metadata['filter'] <- FALSE
  
  # Flag markers with physical location
  geno_usable_idx <- which(!complete.cases(geno_metadata[,c('chrom','pos')]),)
  geno_metadata[geno_usable_idx,'filter'] <- TRUE
  no_pos <- length(geno_usable_idx)
  
  if (no_pos > 0){
    cli::cli_inform("Were found {no_pos} markers
                   without physicall location, flaged")
  }
  
  # Flag colocalized markers
  
  dup_pos_idx <- which(duplicated(geno_metadata[,c('chrom','pos')]))
  no_dup_pos <- length(dup_pos_idx)
  if (no_dup_pos > 0 ){
    cli::cli_inform("Were found {no_dup_pos} markers with duplicated position")
    geno_metadata[dup_pos_idx,'filter'] <- TRUE
  }
  
  # Check duplicated ids
  id_dup = duplicated(geno_metadata['id'])
  
  if(sum(id_dup) > 0){
    cli::cli_inform("Were found {sum(id_dup)} markers with duplicated id. \n
                  Renamed with chrom_pos nomenclature")  
  }
  
  geno_metadata <- geno_metadata %>%
    mutate(
      id = if_else(row_number() %in% id_dup,  paste0(chrom, "_", pos), id),
      allele_count = str_count(alt, ",") + 1,
      ref_len = if_else(str_detect(ref, "^[ACGT]$"),nchar(ref),NA),
      alt_len = if_else(str_detect(alt, "^[ACGT]$"),nchar(str_replace(alt, ",", "")),NA) 
    ) 
  
  # Flag no reference alleles
  ref_alt_na_idx <- which(is.na(geno_metadata['ref_len']) | is.na(geno_metadata['alt_len']))
  no_ref_alt <- length(ref_alt_na_idx)
  
  if(no_ref_alt > 0){
    cli::cli_inform("Were found {no_ref_alt} without ref and alt allele data")
    geno_metadata[ref_alt_na_idx,'filter'] <- TRUE
  }
  
  
    # flag multiallelic 
  multi_allelic_idx <- which(geno_metadata['allele_count'] > 1)
  multi_allelic <- length(multi_allelic_idx)
  
  if(multi_allelic > 0){
    cli::cli_warn("Were found {multi_allelic} multi-allelic markers")
    geno_metadata[multi_allelic_idx,'filter'] <- TRUE
  }

  # flag indels
  
  indel_idx <- which(geno_metadata['ref_len'] > 1)
  indel_idx <- unique(c(indel_idx,which(geno_metadata['alt_len'] > 1)))
  
  
  indel <- length(indel_idx)
  
  if(indel > 0){
    cli::cli_inform("Were found {indel} indel markers")
    geno_metadata[indel_idx,'filter'] <- TRUE
  }
  
  return(geno_metadata)
}

#' From locus genotype call data, get the allelic dosage given the alleles and ploidity
#'
#' This function takes a list of genotype calls, a named vector of allele counts,
#' and the ploidity level as input, and returns a list of allelic dosages for the
#' genotype calls. The allelic dosage is the count of the alternative allele in
#' the genotype call.
#'
#' The function first generates all possible genotype calls for the given alleles
#' and ploidity level using the `get_all_gt_calls` function. It then calculates
#' the allelic dosages for these possible genotype calls using the `convert_gt_to_dosage`
#' function, treating the second allele as the alternative allele.
#'
#' Finally, the function replaces the genotype calls in the input list with their
#' corresponding allelic dosages using the `replace_strings_with_integers` function.
#'
#' If a genotype call in the input list is not found in the set of possible genotype
#' calls, its allelic dosage will be set to NA.
#'
#' @param l A list of genotype calls, e.g., c("AG", "GG", "AA").
#' @param alleles ref and alternative allele, e.g., c("A", "C").
#' @param ploidity Integer. The ploidity level of the organism.
#'
#' @return A list of integers, representing the allelic dosages for the input
#'         genotype calls.
#' @export
#'
#' @examples
#' genotypes <- c("AG", "GG", "AA")
#' allele_counts <- c(A = 10, G = 20)
#' get_allelic_dosage(genotypes, allele_counts, 2)  # Returns list(1, 2, 0)
#'
#' # Example with missing genotype call
#' genotypes <- c("AG", "XX", "AA")
#' allele_counts <- c(A = 10, G = 20)
#' get_allelic_dosage(genotypes, allele_counts, 2)  # Returns list(1, NA, 0)
get_allelic_dosage <- function(l, alleles, ploidity, sep = "") {
  alleles_c <- unlist(stringr::str_split(alleles, "/"))
  # All possible genotype calls
  possible_gt_calls <- get_all_gt_calls(alleles_c, ploidity, sep)
  # Get dosage given the alternative allele
  possible_dosage <- convert_gt_to_dosage(possible_gt_calls, alleles_c[2], ploidity)
  dosages <- replace_strings_with_integers(possible_dosage, l)
  return(dosages)
}


#' Get all possible genotype calls given a unique set of alleles
#'
#' This function generates all possible genotype calls for a given set of alleles
#' and ploidity level. The genotype calls are represented as strings of characters,
#' with each allele being a single character.
#'
#' The function uses a recursive approach to generate all possible combinations of
#' alleles for the specified ploidity level. For example, with two alleles "A" and "B",
#' and a ploidity of 2 (diploid), the function would generate the following genotype
#' calls: "AA", "AB", "BA", "BB".
#'
#' @param alleles List[String]. A list of unique alleles, e.g., c("A", "B", "C").
#' @param ploidity Integer. The ploidity level of the organism.
#'
#' @return List[String]. A list of all possible genotype calls for the given alleles
#'         and ploidity level.
#' @export
#'
#' @examples
#' get_all_gt_calls(c("A", "B"), 2)  # Returns c("AA", "AB", "BA", "BB")
#' get_all_gt_calls(c("A", "B", "C"), 3)  # Returns all 27 possible triploid calls
#' get_all_gt_calls(c("A"), 1)  # Returns c("A")
get_all_gt_calls <- function(alleles, ploidity, sep = "") {
  generate_calls <- function(prefix, ploidity, del = sep) {
    if (ploidity == 0) {
      if(nchar(sep) > 0){
        out <- substr(prefix, 1, nchar(prefix)-1)  
      } else {
        out <- substr(prefix, 1, nchar(prefix))  
      }
      
      return(out)
    }
    calls <- c()
    for (allele in alleles) {
      call <- paste0(prefix, allele,del)
      calls <- c(calls, generate_calls(call, ploidity - 1))
    }
    return(calls)
  }
  generate_calls("", ploidity)
}



#' Given a list of genotype calls, get the dosage of each one
#'
#' This function takes a list of genotype calls, an alternative allele, and the ploidity level
#' of the organism as input, and returns a list of allelic dosages corresponding to each
#' genotype call in the input list.
#'
#' @param locus List. A list of genotype calls, e.g., c("AG", "GG", "AA").
#' @param alt_allele String. The alternative allele, e.g., "A", "G".
#' @param ploidity Integer. The ploidity level of the organism, default is 2 (diploid).
#'
#' @return A list of integers, representing the allelic dosages of the alternative allele
#'         for each genotype call in the input list.
#' @export
#'
#' @examples
#' genotype_calls <- c("AG", "GG", "AA")
#' convert_gt_to_dosage(genotype_calls, "A")  # Returns list(1, 0, 2)
#' convert_gt_to_dosage(genotype_calls, "G")  # Returns list(1, 2, 0)
#' convert_gt_to_dosage(c("AAA", "GGG"), "A", 3)  # Returns list(3, 0) (triploid)
#' convert_gt_to_dosage(c(NA, "AG"), "A")  # Returns list(NA, 1)
convert_gt_to_dosage <- function(locus, alt_allele, ploidity = 2) {
  l <- sapply(locus,
              genocall_to_allelic_dosage,
              alt_allele = alt_allele,
              ploidity = ploidity)
  return(l)
}

#' Genotype call to allelic dosage of alternative allele
#'
#' This function takes a genotype call, an alternative allele, and the ploidity level
#' of the organism as input, and returns the allelic dosage of the alternative allele
#' in the genotype call.
#'
#' The genotype call is expected to be a string of characters representing the alleles,
#' with each allele being a single character. For example, "AG" represents a diploid
#' genotype with one allele being "A" and the other being "G".
#'
#' The allelic dosage is the count of the alternative allele in the genotype call.
#' For example, if the genotype call is "AG" and the alternative allele is "A", the
#' allelic dosage would be 1.
#'
#' If the genotype call is missing (represented as NA or an empty string), the
#' function returns NA.
#'
#' @param genotype_call String. Genotype call, e.g., "AG", "AAA" (for triploid).
#' @param alt_allele String. Alternative allele, e.g., "A", "G".
#' @param ploidity Integer. Ploidity level of the organism, default is 2 (diploid).
#'
#' @return Integer. The allelic dosage of the given genotype call for the alternative allele.
#' @export
#'
#' @examples
#' genocall_to_allelic_dosage("AG", "A")  # Returns 1
#' genocall_to_allelic_dosage("GG", "A")  # Returns 0
#' genocall_to_allelic_dosage("AAA", "A", 3)  # Returns 3 (triploid)
#' genocall_to_allelic_dosage(NA, "A")  # Returns NA
genocall_to_allelic_dosage <- function(genotype_call, alt_allele, ploidity = 2) {
  if (!is.na(nchar(genotype_call))) {
    # Genotype call successfully genotyped
    allele_length <- nchar(genotype_call) / ploidity
    
    # List with each allele as element
    split_genotype <- substring(genotype_call, seq(1, nchar(genotype_call), allele_length),
                                seq(allele_length, nchar(genotype_call), allele_length))
    
    # Matches of alt allele are the dosage
    dosage <- length(which(split_genotype == alt_allele))
    return(dosage)
  } else {
    # Genotype call missed
    return(NA)
  }
}


#' Replace a list of strings with their corresponding integer values
#'
#' This function takes two inputs: a named list or vector with string keys and integer values,
#' and a list of strings to be replaced. It replaces each string in the second list with the
#' corresponding integer value from the first list, based on the string keys.
#'
#' If a string in the second list does not have a corresponding key in the first list,
#' it will be replaced with NA.
#'
#' @param lookup_table A named list or vector with string keys and integer values.
#' @param strings_to_replace A list of strings to be replaced with their corresponding integer values.
#'
#' @return A list of integers, where each string in the input list has been replaced with its
#'         corresponding integer value from the lookup table, or NA if no match was found.
#' @export
#'
#' @examples
#' lookup <- c(A = 1, B = 2, C = 3)
#' strings <- c("B", "A", "D", "C")
#' replace_strings_with_integers(lookup, strings)  # Returns list(2, 1, NA, 3)
replace_strings_with_integers <- function(lookup_table, strings_to_replace) {
  # Use match to find the indices of the strings in the lookup table
  indices <- match(strings_to_replace, names(lookup_table))
  
  # Replace the strings with the corresponding integer values
  # or NA if no match was found
  integer_values <- lookup_table[indices]
  
  return(integer_values)
}


