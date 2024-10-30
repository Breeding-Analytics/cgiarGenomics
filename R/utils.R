#' Read a text table and use common na representations for convert to NA
#'
#' @param path Path to tabular data
#' @param sep  Delimitator of the table cells
#'
#' @return Returns a dataframe
#' @export
#'
#' @examples
read_tabular_geno <- function(path, sep = c('\t',',')){
  
  # Validate params
  sep = match.arg(sep)
  
  # Missing data representation values
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

split_strings <- function(l, ploidity){
  sapply(l, function(x){
    allele_length <- nchar(x) / ploidity
    # List with each allele as element
    split_genotype <- substring(x, seq(1, nchar(x), allele_length),
                                seq(allele_length, nchar(x), allele_length))
    return(split_genotype)
  })
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
get_all_gt_calls <- function(alleles, ploidity) {
  generate_calls <- function(prefix, ploidity) {
    if (ploidity == 0) {
      return(prefix)
    }
    calls <- c()
    for (allele in alleles) {
      call <- paste0(prefix, allele)
      calls <- c(calls, generate_calls(call, ploidity - 1))
    }
    return(calls)
  }
  generate_calls("", ploidity)
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
#' @param allele_count A named vector of allele counts, e.g., c(A = 10, G = 20).
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
get_allelic_dosage <- function(l, allele_count, ploidity) {
  alleles <- names(allele_count)
  
  # All possible genotype calls
  possible_gt_calls <- get_all_gt_calls(alleles, ploidity)
  
  # Get dosage given the alternative allele
  possible_dosage <- convert_gt_to_dosage(possible_gt_calls, alleles[2], ploidity)
  
  dosages <- replace_strings_with_integers(possible_dosage, l)
  return(dosages)
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



