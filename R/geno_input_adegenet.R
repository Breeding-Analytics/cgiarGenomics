

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
read_hapmap <- function(path, ploidity = 2, sep = c("","/","|")) {
  # Validate arguments
  if (!file.exists(path)){
    cli::cli_abort("`path` don't exist. Verify if is writed properly {path}")
  }
  if (!rlang::is_integerish(ploidity)) {
    cli::cli_abort("`ploidity` must be a round number not {ploidity}")
  }
  sep = match.arg(sep)
  
  
  # Read the genotype data from the tabular file
  table <- read_tabular_geno(path)
  
  # Define the expected column names for the hapmap file
  hapmap_snp_attr <- c('rs#', 'alleles', 'chrom', 'pos')
  
  # Check if the first 4 columns in the input file match the expected column names
  if (length(intersect(hapmap_snp_attr, colnames(table)[1:4])) != 4) {
    cli::cli_abort("Hapmap file doesn't have  first four standard column names: {hapmap_snp_attr}, have: {colnames(table)[1:4]}")
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
  
  cli::cli_inform(int_message)
  
  # Check if individual labels and loci IDs are unique
  if (length(unique(individuals)) != length(individuals)) {
    cli::cli_abort("Fatal Error: Individual labels are not unique, check and edit your input file\n")
  }
  
  # Metadata processing, get the first allele as reference and the remaining
  # as alternative
  meta <- table[1:11] %>% 
    mutate(ref = stringr::str_split(alleles, '/', simplify = TRUE)[,1],
           alt = purrr::map_chr(stringr::str_split(alleles, '/'), \(.x) paste(.x[-1], collapse = ","))) %>% 
    mutate(alt = na_if(alt, "")) %>% 
    select('rs#', chrom, pos, ref, alt) %>% 
    rename(id = 'rs#')
  
  meta <- process_metadata(meta)
  
  
  # Read marker columns and transpose to obtain samples x snps
  mt <- t(table[!meta$filter, c(12:dim(table)[2])])
  allele_set <- paste(meta$ref[!meta$filter], meta$alt[!meta$filter], sep='/')
  gt <- mapply(function(col, arg, ploidity) get_allelic_dosage(mt[,col], arg,ploidity),
               col = seq(1,dim(mt)[2]), 
               arg = allele_set,
               ploidity = ploidity)
  
  gl <- new("genlight",
            gt,
            ploidy = ploidity,
            loc.names = meta$id[!meta$filter],
            ind.names = individuals,
            chromosome = meta$chrom[!meta$filter],
            position = meta$pos[!meta$filter])
  adegenet::alleles(gl) <- allele_set
  gl_recalc <- recalc_metrics(gl)
  return(gl)
}

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
read_vcf <- function(path, ploidity = 2, na_reps = c("-", "./.")) {
  
  if (!file.exists(path)){
    cli::cli_abort("`path` don't exist. Verify if is writed properly {path}")
  }
  if (!rlang::is_integerish(ploidity)) {
    cli::cli_abort("`ploidity` must be a round number not {ploidity}")
  }
  # Read the VCF file
  vcf <- vcfR::read.vcfR(path)
  
  # Get the metadata from the VCF file
  meta_vcf <- as.data.frame(vcfR::getFIX(vcf))
  
  meta <- meta_vcf %>% 
    rename(id = ID,
           chrom = CHROM,
           pos = POS,
           ref = REF,
           alt = ALT) %>% 
  select(id, chrom, pos, ref, alt)
  
  meta <- process_metadata(meta)
  mt <- t(vcfR::extract.gt(vcf, return.alleles = T)[!meta$filter,])
  if (length(na_reps) > 0) {
    idx <- which(mt %in% na_reps)
    mt[idx] <- NA
  }
  
  individuals <- rownames(mt)
  
  allele_set <- paste(meta$ref[!meta$filter], meta$alt[!meta$filter], sep='/')
  gt <- mapply(function(col, arg, ploidity, sep) get_allelic_dosage(mt[,col], arg,ploidity,sep),
               col = seq(1,dim(mt)[2]), 
               arg = allele_set,
               ploidity = ploidity,
               sep = "/")
  
  gl <- new("genlight",
            gt,
            ploidy = ploidity,
            loc.names = meta$id[!meta$filter],
            ind.names = individuals,
            chromosome = meta$chrom[!meta$filter],
            position = meta$pos[!meta$filter])
  adegenet::alleles(gl) <- allele_set
  gl_recalc <- recalc_metrics(gl)
  return(gl_recalc)
}


#' read_DArTSeq_SNP
#'
#' This function reads a DArTSeq SNP dataset from a CSV file and converts it into a genlight object using
#' the DartR library reading functions.
#'
#' @param path Path to the DArTSeq CSV file.
#' @param snp_id Column name of the SNP ID data in the CSV file.
#' @param chr_name Column name of the chromosome where the SNP is located in the CSV file.
#' @param pos_name Column name of the physical position where the SNP is located in the CSV file.
#'
#' @return A genlight object with mapped SNPs.
#' @export
#'
#' @examples
#' read_DArTSeq_SNP("path/to/dartseq/file.csv", snp_id = "SnpID", chr_name = "Chr", pos_name = "Position")
read_DArTSeq_SNP <- function(path, snp_id, chr_name, pos_name) {
  if (!file.exists(path)){
    cli::cli_abort("`path` don't exist. Verify if is writed properly {path}")
  }
  
  # Read the DArTSeq CSV file
  gl <- dartR::gl.read.dart(path)
  
  
  
  
  if(!hasName(gl@other$loc.metrics, snp_id)){
    cli::cli_abort("The input snp_id: {snp_id} column doesn't exist in the dartSeq file")
  }
  
  if(!hasName(gl@other$loc.metrics, chr_name)){
    cli::cli_abort("The input chr_name: {chr_name} column doesn't exist in the dartSeq file")
  }
  
  if(!hasName(gl@other$loc.metrics, pos_name)){
    cli::cli_abort("The input pos_name: {pos_name} column doesn't exist in the dartSeq file")
  }
  
  # DArTSeq stores position and chromosome information in the 'other' slot
  # It is necessary to modify it to be compatible with the native genlight object
  pos <- c(gl@other$loc.metrics[, pos_name])
  chrom <- c(gl@other$loc.metrics[, chr_name])
  
  # Remove SNPs where the tag is unmapped (no position data)
  no_pos <- which(is.na(pos) | pos <= 0)
  no_chrom <- which(is.na(chrom))
  no_loc <- unique(c(no_pos, no_chrom))
  
  if(length(no_loc) > 0){
    cli::cli_inform("Were detected {length(no_loc)} \n
                    SNPs withoud position data, removed.")
  }
  
  # Remove SNPs without position data
  gl <- gl[, -c(no_loc)]
  
  # Assign the position and chromosome information to the genlight object
  adegenet::position(gl) <- gl@other$loc.metrics[, pos_name]
  adegenet::chromosome(gl) <- gl@other$loc.metrics[, chr_name]
  
  gl_recalc <- recalc_metrics(gl)
  return(gl_recalc)
  
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
#'
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


