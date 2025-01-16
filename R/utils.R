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
  
  # Filtering modifications
  filt_mods <- modifications %>% 
    filter(!grepl("^imputation", reason))
  
  
  ind_out <- filt_mods %>% 
    dplyr::filter(is.na(col)) %>% 
    dplyr::pull(row)
  
  loc_out <- filt_mods %>% 
    dplyr::filter(is.na(row)) %>% 
    dplyr::pull(col)
  
  ind_idx <- seq(1:length(adegenet::indNames(gl)))[-ind_out]
  loc_idx <- seq(1:length(adegenet::locNames(gl)))[-loc_out]
  

  filt_gl <- gl[ind_idx, loc_idx]
  return(filt_gl)
  
  imp_mods <- modifications %>% 
    filter(grepl("^imputation", reason))
  
  mt <- as.matrix(filt_gl)

  mt[cbind(imp_mods$row, imp_mods$col)] <- imp_mods$value
  

  imp_gl <- new("genlight",
                mt,
                loc.names = filt_gl@loc.names,
                ind.names = filt_gl@ind.names,
                chromosome = filt_gl@chromosome,
                position = filt_gl@position)
  
  adegenet::alleles(imp_gl) <- adegenet::alleles(filt_gl)
  imp_gl <- recalc_metrics(imp_gl)
  
  return(imp_gl)
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





