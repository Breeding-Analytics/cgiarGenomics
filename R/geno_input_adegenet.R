

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

  gt <- mapply(function(col, arg, ploidity, sep) get_allelic_dosage(mt[,col], arg,ploidity, sep),
               col = seq(1,dim(mt)[2]), 
               arg = allele_set,
               ploidity = ploidity,
               sep = sep)
  
  gl <- new("genlight",
            gt,
            ploidy = ploidity,
            loc.names = meta$id[!meta$filter],
            ind.names = individuals,
            chromosome = meta$chrom[!meta$filter],
            position = meta$pos[!meta$filter])
  adegenet::alleles(gl) <- allele_set
  gl <- recalc_metrics(gl)
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
read_vcf <- function(path, ploidity = 2, na_reps = c("-", "./."), sep="/") {
  
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
  
  allele_set <- paste(meta$ref[!meta$filter], meta$alt[!meta$filter], sep=sep)

  gt <- mapply(function(col, arg, ploidity, sep) get_allelic_dosage(mt[,col], arg,ploidity, sep),
               col = seq(1,dim(mt)[2]), 
               arg = allele_set,
               ploidity = ploidity,
               sep = sep)
  
  gl <- new("genlight",
            gt,
            ploidy = ploidity,
            loc.names = meta$id[!meta$filter],
            ind.names = individuals,
            chromosome = meta$chrom[!meta$filter],
            position = meta$pos[!meta$filter])
  adegenet::alleles(gl) <- allele_set
  gl <- recalc_metrics(gl)
  return(gl)
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
  
  gl <- recalc_metrics(gl)
  return(gl)
  
}

#' extract_dart_gt
#'
#' @param data Dataframe. Dataframe with Allelic dosage data
#' @param sample.name.row Integer. Row number which the sample names are (default = 2)
#' @param first.data.row Integer. Row number of the first genotype call
#' @param first.data.col Integer. Col number of the first genotype call
#'
#' @return Dataframe. Returns a Dataframe with the allelic dosage data converted
#' to integers and the column names as the samples.
#' @export
#'
#' @examples
extract_dart_gt <- function(data,
                            sample.name.row= 2,
                            first.data.row = 9,
                            first.data.col=6){
  
  if (!is.data.frame(data)) {
    cli::cli_abort("`sample.name.row` must be a round number not {sample.name.row}")
  }
  if (!rlang::is_integerish(sample.name.row)) {
    cli::cli_abort("`sample.name.row` must be a round number not {sample.name.row}")
  }
  if (!rlang::is_integerish(first.data.row)) {
    cli::cli_abort("`first.data.row` must be a round number not {first.data.row}")
  }
  if (!rlang::is_integerish(first.data.col)) {
    cli::cli_abort("`first.data.col` must be a round number not {first.data.col}")
  }
  
  # rows with data
  rows <- first.data.row:nrow(data)
  cols <- first.data.col:ncol(data)
  # sample names
  id <- as.character(data[sample.name.row,cols])
  # Samples with duplicated id?
  dupes <- unique(id[which(duplicated(id))])
  if(length(dupes)>0){
    cli::cli_abort("The samples: {dupes} are duplicated, fix it.")
  }
  out_data <- data[rows, cols]
  colnames(out_data) <- id
  return(out_data)
}

#' read_one_row_dart_gt_data
#' 
#' This function takes as input a path a DArT file where each variant is coded
#' in a single line. It splits the metadata and the allelic dosages.
#'
#' @param path 
#' @param sample.name.row 
#' @param first.data.row 
#' @param first.data.col 
#' @param snp_id_col 
#'
#' @return
#' @export
#'
#' @examples
read_one_row_dart_gt_data <- function(path,
                                      sample.name.row= 2,
                                      first.data.row = 9,
                                      first.data.col=6,
                                      snp_id_col){
  if (!file.exists(path)){
    cli::cli_abort("`path` don't exist. Verify if is writed properly {path}")
  }
  
  data_raw <- read.csv(path, header=F, check.names=F)
  data <- extract_dart_gt(data_raw,
                          sample.name.row,
                          first.data.row,
                          first.data.col)
  
  
  # Markers with duplicated id
  meta_data <- data_raw[first.data.row:nrow(data_raw), 1:first.data.col-1]
  colnames(meta_data) <- data_raw[first.data.row-1, 1:first.data.col-1]
  
  if(!hasName(meta_data, snp_id_col)){
    cli::cli_abort("The input `snp_id_col`: {snp_id_col} column doesn't exist in the file")
  }
  
  marker_ids <- as.vector(unlist(meta_data[snp_id_col]))
  l_dupes <- unique(marker_ids[which(duplicated(marker_ids))])
  
  geno <- apply(data,2,as.integer)
  
  geno <- as.data.frame(geno)
  geno$marker_id <- marker_ids
  
  return(list(gt = geno, dosage_meta = meta_data))
}

parse_counts_metadata <- function(path,
                                  sample.name.row= 7,
                                  first.data.row = 9,
                                  first.data.col=6,
                                  snp_id_col = 'CloneID',
                                  allele_id_col = 'AlleleID',
                                  allele_sep = "\\|",
                                  allele_sequence_col = 'AlleleSequence',
                                  ref_name = 'Ref',
                                  alt_name = 'Alt'){
  
  data_raw <- read.csv(path, header=F, check.names=F)
  
  meta_data <- data_raw[first.data.row:nrow(data_raw), 1:first.data.col-1]
  colnames(meta_data) <- data_raw[first.data.row-1, 1:first.data.col-1]
  
  if(!hasName(meta_data, snp_id_col)){
    cli::cli_abort("The input `snp_id_col`: {snp_id_col} column doesn't exist in the file")
  }
  if(!hasName(meta_data, allele_id_col)){
    cli::cli_abort("The input `allele_id_col`: {allele_id_col} column doesn't exist in the file")
  }
  if(!hasName(meta_data, allele_sequence_col)){
    cli::cli_abort("The input `allele_sequence_col`: {allele_sequence_col} column doesn't exist in the file")
  }
  
  meta_data <- meta_data %>% 
    mutate(allele_type = stringr::str_split_fixed(!!sym(allele_id_col), allele_sep, 2)[, 2])
  
  # Wider the metadata table
  meta_data_wide <- meta_data %>% 
    select(!!sym(snp_id_col),
           !!sym(allele_sequence_col),
           allele_type) %>% 
    tidyr::pivot_wider(
      id_cols = !!sym(snp_id_col),
      names_from = allele_type,
      values_from = !!sym(allele_sequence_col)
    )
  
  alleles <- purrr::map2_dfr(as.vector(unlist(meta_data_wide[,ref_name])),
                             as.vector(unlist(meta_data_wide[,alt_name])),
                             \(.x,.y) compare_tags(.x, .y))
  
  wide_alleles <- cbind(meta_data_wide, alleles)
  
  
  return(wide_alleles)
}

read_DArTag_count_dosage <- function(dosage_path,
                                     counts_path,
                                     ploidity,
                                     sample.name.row = 7,
                                     first.data.row = 9,
                                     first.data.col = 6,
                                     snp_id_col_dosage = 'MarkerID',
                                     snp_id_col_count = 'CloneID',
                                     allele_id_col = 'AlleleID',
                                     allele_sep = "\\|",
                                     allele_sequence_col = 'AlleleSequence',
                                     ref_name = 'Ref',
                                     alt_name = 'Alt',
                                     dosage_col_chr = 'Chrom',
                                     dosage_col_pos = 'ChromPos'){
  
  if (!file.exists(dosage_path)){
    cli::cli_abort("`dosage_path` don't exist. Verify if is writed properly {dosage_path}")
  }
  if (!file.exists(counts_path)){
    cli::cli_abort("`counts_path` don't exist. Verify if is writed properly {counts_path}")
  }
  # Parse metadata (Allele sequences)
  meta_data <- parse_counts_metadata(
    counts_path,
    sample.name.row,
    first.data.row,
    first.data.col,
    snp_id_col_count,
    allele_id_col,
    allele_sep,
    allele_sequence_col,
    ref_name,
    alt_name)
  
  # Parse genotype call information
  dosage_data <- read_one_row_dart_gt_data(
    dosage_path, 
    sample.name.row,
    first.data.row,
    first.data.col,
    snp_id_col_dosage)
  
  gt <- dosage_data$gt
  
  dosage_meta <- dosage_data$dosage_meta
  
  if(!hasName(dosage_meta, dosage_col_chr)){
    cli::cli_abort("The input `dosage_col_chr`: {dosage_col_chr} column doesn't exist in the file")
  }
  if(!hasName(dosage_meta, dosage_col_pos)){
    cli::cli_abort("The input `dosage_col_pos`: {dosage_col_pos} column doesn't exist in the file")
  }
  
  
  dosage_meta[,dosage_col_pos] <- as.integer(dosage_meta[,dosage_col_pos])
  
  # Check if markers reported in dosage are in counts 
  
  orphan_snps <- which(!gt$marker_id %in% meta_data[, snp_id_col_count])
  if(length(orphan_snps) > 0){
    cli::cli_warn("The markers: {gt$marker_id[orphan_snps]} without Allele sequence data, removed.")
  }
  
  
  merged_data <- merge(gt, meta_data,
                       by.x = "marker_id",
                       by.y = snp_id_col_count)
  
  
  
  merged_data <- merge(merged_data, dosage_meta,
                       by.x = "marker_id",
                       by.y = snp_id_col_dosage)
  
  
  nrows <- dim(gt)[1]
  ncols <- dim(gt)[2]
  individuals <- colnames(gt)[which(colnames(gt) != "marker_id")]
  
  
  loci <- paste0(c(merged_data[1:5,dosage_col_chr]), "_", c(merged_data[1:5,dosage_col_pos]))
  
  report_bullets <- cli::cli_bullets(c("Input data should have loci as rows and individuals as columns.",
                                       "i" = "{nrows} Loci, first 5: {loci}",
                                       "i" = "{length(individuals)} Individuals, first 5: {individuals[1:5]}"))
  
  cli::cli_inform(report_bullets)
  
  # Metadata processing, get the first allele as reference and the remaining
  # as alternative
  meta <- merged_data[,c('marker_id', dosage_col_chr, dosage_col_pos, 'ref_allele', 'alt_allele')]
  colnames(meta) <- c('id', 'chrom', 'pos', 'ref', 'alt')
  
  
  meta <- process_metadata(meta)
  
  # Read marker columns and transpose to obtain samples x snps
  mt <- t(merged_data[!meta$filter, individuals])
  allele_set <- paste(meta$ref[!meta$filter], meta$alt[!meta$filter], sep='/')
  
  gl <- new("genlight",
            mt,
            ploidy = ploidity,
            loc.names = meta$id[!meta$filter],
            ind.names = individuals,
            chromosome = meta$chrom[!meta$filter],
            position = meta$pos[!meta$filter])
  adegenet::alleles(gl) <- allele_set
  gl <- recalc_metrics(gl)
  return(gl)
}

compare_tags <- function(ref, alt){
  bases <- c("A", "C", "G", "T")
  hap.ref <- unlist(strsplit(ref, split = ""))
  hap.alt <- unlist(strsplit(alt, split = ""))
  k <- which(hap.ref != hap.alt & hap.ref %in% bases & 
               hap.alt %in% bases)
  if (length(k) == 1) {
    return(c(ref_allele = hap.ref[k], alt_allele = hap.alt[k]))
  }
  else {
    return(c(ref_allele = "N", alt_allele = "N"))
  }
}


