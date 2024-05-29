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
  # Convert genotypes to 0, 1, 2, and NA
  for (i in 1:dim(mt)[2]) {
    v1 <- mt[, i]
    allele_count = get_alleles_count_char(v1)
    if (length(names(allele_count)) >= 3) {
      non_biallelic <- c(i, non_biallelic)
    } else {
      homRef <- paste0(names(allele_count)[1], names(allele_count)[1])
      homAlt <- paste0(names(allele_count)[2], names(allele_count)[2])
      het1 <- paste0(names(allele_count)[1], names(allele_count)[2])
      het2 <- paste0(names(allele_count)[2], names(allele_count)[1])
      mt[, i] <- gsub(homRef, "0", mt[, i])
      mt[, i] <- gsub(homAlt, "2", mt[, i])
      mt[, i] <- gsub(het1, "1", mt[, i])
      mt[, i] <- gsub(het2, "1", mt[, i])
    }
  }
  
  if (length(non_biallelic) > 0) {
    bi_message <- paste0("From ", nrows, "loci, ",
                         length(non_biallelic), " aren't bi-allelic")
    gl <- new("genlight",
              mt[, -c(non_biallelic)],
              ploidy = ploidity,
              loc.names = loci[-c(non_biallelic)],
              ind.names = individuals,
              chromosome = table[-c(non_biallelic), 3],
              position = table[-c(non_biallelic), 4])
  } else {
    bi_message <- "Data confirmed bi-allelic"
    gl <- new("genlight",
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

# Test
# ploidity = 2
# sep = ""
# path <- 'tests/hmp_fmt/geno_pruned_hmp.txt'
# geno <- read_hapmap(path = path)


#' read_vcf
#'
#' This function reads a VCF file (compressed or uncompressed) and converts it into a genlight object.
#'
#' @param path String. Path to the VCF file. It could be compressed.
#' @param ploidity Integer. Ploidity level of the organism. (Default = 2)
#' @param na_reps Vector. A vector containing the NA representations of genotype calls (default: empty).
#'
#' @return A genlight object.
#' @export
#'
#' @examples
#' read_vcf("path/to/vcf/file.vcf", ploidity = 2)
#' read_vcf("path/to/vcf/file.vcf.gz", ploidity = 4, na_reps = c("-", "./."))
read_vcf <- function(path, ploidity = 2, na_reps = c()) {
  # Read the VCF file
  vcf <- vcfR::read.vcfR(path)
  
  # Get the metadata from the VCF file
  meta_vcf <- vcfR::getFIX(vcf)
  
  # Extract the genotype data and transpose it (samples x snps)
  mt <- t(vcfR::extract.gt(vcf, return.alleles = T))
  
  # If there are any NA representations provided, replace them with NA
  if (length(na_reps) > 0) {
    idx <- which(mt %in% na_reps)
    mt[idx] <- NA
  }
  
  # Remove markers where all samples have missing data
  cols_to_remove <- colSums(is.na(mt)) == nrow(mt)
  na_markers <- c(which(cols_to_remove))
  
  if (length(na_markers) > 0) {
    na_message <- paste0("From ", dim(mt)[2], "loci, ",
                         length(na_markers), " have complete missing data, removed.")
    print_log_message(na_message)
    mt <- mt[, !cols_to_remove]
  }
  
  # Create loci IDs from chromosome and position
  loci <- paste0(meta_vcf[!cols_to_remove, "CHROM"], "_", meta_vcf[!cols_to_remove, "POS"])
  individuals <- rownames(mt)
  non_biallelic <- c()
  
  # Convert genotypes to 0, 1, 2, and NA
  for (i in 1:dim(mt)[2]) {
    v1 <- mt[, i]
    allele_count <- get_alleles_count_char(v1, ploidity = ploidity, sep = '/')
    
    # If the marker is not bi-allelic, add it to the non_biallelic vector
    if (length(names(allele_count)) > 3) {
      print(allele_count)
      non_biallelic <- c(i, non_biallelic)
    } else {
      # Substitute genotype codes with 0, 1, 2 for homozygous reference, heterozygous, and homozygous alternative
      homRef <- paste0(names(allele_count)[1], names(allele_count)[1])
      homAlt <- paste0(names(allele_count)[2], names(allele_count)[2])
      het1 <- paste0(names(allele_count)[1], names(allele_count)[2])
      het2 <- paste0(names(allele_count)[2], names(allele_count)[1])
      mt[, i] <- gsub(homRef, 0, mt[, i])
      mt[, i] <- gsub(homAlt, 0, mt[, i])
      mt[, i] <- gsub(het1, 1, mt[, i])
      mt[, i] <- gsub(het2, 1, mt[, i])
    }
  }
  
  if (length(non_biallelic) > 0) {
    bi_message <- paste0("From ", dim(mt)[2], "loci, ",
                         length(non_biallelic), " aren't bi-allelic")
    
    # Remove non-biallelic and markers with missing data
    non_biallelic <- c(non_biallelic, na_markers)
    
    # Create the genlight object with the remaining markers
    gl <- new("genlight",
              mt[, -c(non_biallelic)],
              ploidy = ploidity,
              loc.names = loci[-c(non_biallelic)],
              ind.names = individuals,
              chromosome = meta_vcf[-c(non_biallelic), "CHROM"],
              position = as.numeric(meta_vcf[-c(non_biallelic), "POS"]))
  } else {
    bi_message <- "Data confirmed bi-allelic"
    
    # Create the genlight object with all markers
    gl <- new("genlight",
              mt,
              ploidy = ploidity,
              loc.names = loci,
              ind.names = individuals,
              chromosome = meta_vcf[-c(na_markers), "CHROM"],
              position = meta_vcf[-c(na_markers), "POS"])
  }
  
  print_log_message(bi_message)
  return(gl)
}

#' read_DArTSeq_SNP
#'
#' This function reads a DArTSeq SNP dataset from a CSV file and converts it into a genlight object using
#' the DartR library reading functions.
#'
#' @param dart_path Path to the DArTSeq CSV file.
#' @param snp_id Column name of the SNP ID data in the CSV file.
#' @param chr_name Column name of the chromosome where the SNP is located in the CSV file.
#' @param pos_name Column name of the physical position where the SNP is located in the CSV file.
#'
#' @return A genlight object with mapped SNPs.
#' @export
#'
#' @examples
#' read_DArTSeq_SNP("path/to/dartseq/file.csv", snp_id = "SnpID", chr_name = "Chr", pos_name = "Position")
read_DArTSeq_SNP <- function(dart_path, snp_id, chr_name, pos_name) {
  # Read the DArTSeq CSV file
  gl <- dartR::gl.read.dart(dart_path)
  
  # DArTSeq stores position and chromosome information in the 'other' slot
  # It is necessary to modify it to be compatible with the native genlight object
  pos <- c(gl@other$loc.metrics[, pos_name])
  chrom <- c(gl@other$loc.metrics[, chr_name])
  
  # Remove SNPs where the tag is unmapped (no position data)
  no_pos <- which(is.na(pos) | pos <= 0)
  no_chrom <- which(is.na(chrom))
  no_loc <- unique(c(no_pos, no_chrom))
  
  message <- paste("Were detected:",
                   length(no_loc),
                   "SNPs without position data, removed.")
  print_log_message(message)
  
  # Remove SNPs without position data
  gl <- gl[, -c(no_loc)]
  
  # Assign the position and chromosome information to the genlight object
  position(gl) <- gl@other$loc.metrics[, pos_name]
  chromosome(gl) <- gl@other$loc.metrics[, chr_name]
  
  return(gl)
}

#' Read a DArTSeq Presence/Absence file
#'
#' This function reads a DArTSeq Presence/Absence (PA) dataset from a CSV file and converts it into a genlight object using
#' the DartR library reading functions.
#'
#' @param dart_path String, Path to the DArTSeq PA CSV file.
#' @param marker_id String, Column name in the CSV file that represents the marker ID.
#' @param chr_name String, Column name in the CSV file that represents the chromosome name.
#' @param pos_name String, Column name in the CSV file that represents the physical position.
#'
#' @return A genlight object.
#' @export
#'
#' @examples
#' read_DArTSeq_PA("path/to/dartseq/pa/file.csv", marker_id = "MarkerID", chr_name = "Chr", pos_name = "Position")
read_DArTSeq_PA <- function(dart_path, marker_id, chr_name, pos_name) {
  # Read the DArTSeq PA CSV file
  gl <- dartR::gl.read.silicodart(dart_path)
  
  # DArTSeq stores position and chromosome information in the 'other' slot
  # It is necessary to modify it to be compatible with the native genlight object
  pos <- c(gl@other$loc.metrics[, pos_name])
  chrom <- c(gl@other$loc.metrics[, chr_name])
  
  # Remove markers where the tag is unmapped (no position data)
  no_pos <- which(is.na(pos) | pos <= 0)
  no_chrom <- which(is.na(chrom))
  no_loc <- unique(c(no_pos, no_chrom))
  
  message <- paste("Were detected:",
                   length(no_loc),
                   "SNPs without position data, removed.")
  print_log_message(message)
  
  # Remove markers without position data
  gl <- gl[, -c(no_loc)]
  
  # Assign the position and chromosome information to the genlight object
  position(gl) <- gl@other$loc.metrics[-c(no_loc), pos_name]
  chromosome(gl) <- gl@other$loc.metrics[-c(no_loc), chr_name]
  
  return(gl)
}

#' Read a DArTTag file
#'
#' This function reads a DArTTag report (counts and dosage CSV files) and converts it into a genlight object.
#' It first converts the DArTTag report to a VCF file using the polyBreedR library, and then reads the VCF file
#' into a genlight object using the read_vcf function.
#'
#' @param counts.file String. Path to the DArTTag counts CSV file.
#' @param dosage.file String. Path to the DArTTag dosage CSV file.
#' @param ploidity Integer, ploidity level of the organism (default: 4).
#' @param na_reps Vector. A vector containing the NA representations of genotype calls (default: ".").
#'
#' @return A genlight object.
#' @export
#'s
#' @examples
#' read_DArT_Tag("path/to/counts.csv", "path/to/dosage.csv")
#' read_DArT_Tag("path/to/counts.csv", "path/to/dosage.csv", ploidity = 2, na_reps = c("-", "./."))
read_DArT_Tag <- function(counts.file, dosage.file, ploidity = 4, na_reps = c(".")) {
  # Convert the DArTTag report to a VCF file
  polyBreedR::dart2vcf(counts.file = counts.file,
                       dosage.file = dosage.file,
                       ploidy = ploidity,
                       vcf.file = "tmp.vcf.gz")
  
  # Read the VCF file into a genlight object
  gl <- read_vcf("tmp.vcf.gz", ploidity = ploidity, na_reps = na_reps)
  
  return(gl)
}
