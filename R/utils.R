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

#' Write genlight object to VCF file
#'
#' This function converts a genlight object to a VCF file and writes it to an output path.
#'
#' @param gl A genlight object.
#' @param output_path String. Path to the output VCF file. It could end in `.gz` for compression.
#' @param na_rep String. NA representation for genotypes in the VCF file. Default is "./.".
#'
#' @export
write_vcf <- function(gl, output_path, na_rep = "./.") {
  if (!inherits(gl, "genlight")) {
    cli::cli_abort("`gl` must be a genlight object.")
  }
  if (!is.character(output_path) || length(output_path) != 1) {
    cli::cli_abort("`output_path` must be a single character string.")
  }
  
  n_loci <- length(adegenet::locNames(gl))
  n_samples <- length(adegenet::indNames(gl))
  
  gt_matrix <- as.matrix(gl)
  
  # Extract alleles
  alleles <- adegenet::alleles(gl)
  if (is.null(alleles)) {
    ref_alleles <- rep("N", n_loci)
    alt_alleles <- rep("N", n_loci)
  } else {
    split_alleles <- strsplit(as.character(alleles), "/")
    ref_alleles <- sapply(split_alleles, function(x) {
      if (length(x) >= 1 && !is.na(x[1])) x[1] else "N"
    })
    alt_alleles <- sapply(split_alleles, function(x) {
      if (length(x) >= 2 && !is.na(x[2])) x[2] else "N"
    })
  }
  
  # Ploidy
  p_ind <- adegenet::ploidy(gl)
  if (length(p_ind) == 1) {
    p_ind <- rep(p_ind, n_samples)
  }
  
  # Translate dosages to VCF genotype format
  V <- matrix("", nrow = n_samples, ncol = n_loci)
  for (i in seq_len(n_samples)) {
    p <- p_ind[i]
    possible_dosages <- 0:p
    gt_strings <- sapply(possible_dosages, function(d) {
      paste(c(rep("0", p - d), rep("1", d)), collapse = "/")
    })
    names(gt_strings) <- as.character(possible_dosages)
    
    row_na_rep <- if (is.null(na_rep)) {
      paste(rep(".", p), collapse = "/")
    } else {
      na_rep
    }
    
    row_vals <- gt_matrix[i, ]
    row_chars <- as.character(row_vals)
    translated <- gt_strings[row_chars]
    translated[is.na(row_vals)] <- row_na_rep
    V[i, ] <- translated
  }
  
  # Calculate AC and AN
  ac_vals <- colSums(gt_matrix, na.rm = TRUE)
  if (length(unique(p_ind)) == 1) {
    an_vals <- colSums(!is.na(gt_matrix)) * p_ind[1]
  } else {
    an_vals <- colSums(!is.na(gt_matrix) * p_ind)
  }
  info_val <- paste0("AC=", ac_vals, ";AN=", an_vals)
  
  # Chromosome and position
  chrom_val <- if (!is.null(adegenet::chromosome(gl))) as.character(adegenet::chromosome(gl)) else rep(".", n_loci)
  pos_val <- if (!is.null(adegenet::position(gl))) as.character(adegenet::position(gl)) else rep(".", n_loci)
  
  chrom_val[is.na(chrom_val)] <- "."
  pos_val[is.na(pos_val)] <- "."
  
  id_val <- adegenet::locNames(gl)
  if (is.null(id_val)) id_val <- paste0("Locus_", seq_len(n_loci))
  
  # VCF data columns
  vcf_data <- data.frame(
    `#CHROM` = chrom_val,
    POS = pos_val,
    ID = id_val,
    REF = ref_alleles,
    ALT = alt_alleles,
    QUAL = rep(".", n_loci),
    FILTER = rep("PASS", n_loci),
    INFO = info_val,
    FORMAT = rep("GT", n_loci),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Genotype columns
  gt_cols <- t(V)
  colnames(gt_cols) <- adegenet::indNames(gl)
  vcf_df <- cbind(vcf_data, gt_cols)
  
  # Headers
  header <- c(
    "##fileformat=VCFv4.0",
    "##FILTER=<ID=PASS,Description=\"All filters passed\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">",
    "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">"
  )
  
  # Write
  con <- if (grepl("\\.gz$", output_path, ignore.case = TRUE)) {
    gzfile(output_path, "wt")
  } else {
    file(output_path, "wt")
  }
  writeLines(header, con = con)
  write.table(vcf_df, file = con, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  close(con)
}

check_method_args <- function(dots, allowed, required, method) {
  supplied <- names(dots)
  
  if (is.null(supplied)) {
    supplied <- character()
  }
  
  unnamed <- supplied == ""
  
  if (any(unnamed)) {
    stop(
      "All method-specific arguments passed through '...' must be named.",
      call. = FALSE
    )
  }
  
  unknown <- setdiff(supplied, allowed)
  
  if (length(unknown) > 0) {
    stop(
      sprintf(
        "Unknown argument(s) for method '%s': %s. Allowed arguments are: %s.",
        method,
        paste(unknown, collapse = ", "),
        paste(allowed, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  missing_required <- setdiff(required, supplied)
  
  if (length(missing_required) > 0) {
    stop(
      sprintf(
        "Missing required argument(s) for method '%s': %s.",
        method,
        paste(missing_required, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}
