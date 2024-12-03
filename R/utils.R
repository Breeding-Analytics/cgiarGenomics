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



apply_bioflow_modifications <- function(gl, modifications){
  
  ind_out <- modifications %>% 
    dplyr::filter(margin == 'ind') %>% 
    dplyr::pull(filt_out)
  
  loc_out <- modifications %>% 
    dplyr::filter(margin == 'loc') %>% 
    dplyr::pull(filt_out)
  
  ind_idx <- which(!ind_out %in% adegenet::indNames(gl))
  loc_idx <- which(!loc_out %in% adegenet::locNames(gl))
  
  print(ind_idx)
  print(loc_idx)
  return(gl[ind_idx, loc_idx])
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





