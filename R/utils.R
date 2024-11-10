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

split_strings <- function(l, ploidity){
  sapply(l, function(x){
    allele_length <- nchar(x) / ploidity
    # List with each allele as element
    split_genotype <- substring(x, seq(1, nchar(x), allele_length),
                                seq(allele_length, nchar(x), allele_length))
    return(split_genotype)
  })
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



