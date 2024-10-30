# path <- "tests/DartSeq_fmt/DartSeq/Report_DCob23-8528_SNP_mapping_2.csv"
# snp_id <- "AlleleID"
# pos <- "ChromPosSnp_Common_bean_Phaseolus_acutifolius_v10"
# chrom <- "Chrom_Common_bean_Phaseolus_acutifolius_v10"
# gl <- read_DArTSeq_SNP(path, snp_id, chrom, pos)
# 
# path <- "tests/DartSeq_fmt/DartSeq/Report_DCob23-8528_SilicoDArT_1.csv"
# snp_id <- "CloneID"
# chrom <- "Chrom_Common_bean_Phaseolus_acutifolius_v10"
# pos <- "ChromPosTag_Common_bean_Phaseolus_acutifolius_v10"
# 
# gl2 <- read_DArTSeq_PA(path, snp_id, chrom, pos)
# 



#' read_hapmap
#'
#' This function reads a hapmap file (uncompressed) and converts it into a genlight object.
#'
#' @param path Path to the hapmap file (uncompressed).
#' @param ploidity Integer, ploidity level of the organism (default: 2).
#' @param sep String, separator that divides the alleles (default: "").
#'
#' @return A genlight object.
#' @export
#'
#' @examples
#' read_hapmap("path/to/hapmap/file.txt")
#' read_hapmap("path/to/hapmap/file.txt", ploidity = 4)
#' read_hapmap("path/to/hapmap/file.txt", sep = "|")
read_hapmap <- function(path, ploidity = 2, sep = "") {
  # Read the genotype data from the tabular file
  table <- read_tabular_geno(path)
  
  # Define the expected column names for the hapmap file
  hapmap_snp_attr <- c('rs#', 'alleles', 'chrom', 'pos', 'strand', 'assembly#',
                       'center', 'protLSID', 'assayLSID', 'panel', 'QCcode')
  
  # Check if the first 11 columns in the input file match the expected column names
  if (length(intersect(hapmap_snp_attr, colnames(table)[1:11])) != 11) {
    print_log_message("Hapmap file doesn't have the standard column names")
    stop()
  }
  
  nrows <- dim(table)[1]
  ncols <- dim(table)[2]
  individuals <- colnames(table)[12:ncols]
  loci <- paste0(c(table[1:nrows, 3]), "_", c(table[1:nrows, 4]))
  
  # Print a data integrity message
  int_message <- paste0(
    "Input data should have loci as rows and individuals as columns \n",
    nrows - 1, " Loci, confirming first 5: \n",
    paste(loci[1:5], collapse = " "), "\n",
    ncols - length(hapmap_snp_attr), " Individuals, confirming first 5:\n",
    paste(individuals[1:5], collapse = " ")
  )
  print_log_message(int_message)
  
  # Check if individual labels and loci IDs are unique
  if (length(unique(individuals)) != length(individuals)) {
    cat(error("Fatal Error: Individual labels are not unique, check and edit your input file\n"))
    stop()
  }
  if (length(unique(loci)) != length(loci)) {
    cat(error("Fatal Error: loci not unique, check and edit your input file\n"))
    stop()
  }
  
  # Read marker columns and transpose to obtain samples x snps
  mt <- t(table[, c(12:dim(table)[2])])
  non_biallelic <- c()
  for (i in 1:dim(mt)[2]) {
    v1 <- mt[, i]
    allele_count = get_alleles_count_char(v1, ploidity)
    if (length(names(allele_count)) >= 3) {
      non_biallelic <- c(i, non_biallelic)
    } else {
      alt_allele <- names(allele_count)[2]
      l <- get_allelic_dosage(mt[,i], allele_count, ploidity)
      mt[,i] <- get_allelic_dosage(mt[,i], allele_count, ploidity)
      
    }
  }
  
  if (length(non_biallelic) > 0) {
    bi_message <- paste0("From ", nrows, "loci, ",
                         length(non_biallelic), " aren't bi-allelic")
    gl <- new("genlight",
              package = 'adegenet',
              mt[, -c(non_biallelic)],
              ploidy = ploidity,
              loc.names = loci[-c(non_biallelic)],
              ind.names = individuals,
              chromosome = table[-c(non_biallelic), 3],
              position = table[-c(non_biallelic), 4])
  } else {
    bi_message <- "Data confirmed bi-allelic"
    gl <- new("genlight",
              package = 'adegenet',
              mt,
              ploidy = ploidity,
              loc.names = loci,
              ind.names = individuals,
              chromosome = table[1:nrows, 3],
              position = table[1:nrows, 4])
  }
  print_log_message(bi_message)
  return(gl)
}