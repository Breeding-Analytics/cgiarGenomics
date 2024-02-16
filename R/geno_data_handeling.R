
# Read functions ----------------------------------------------------------

read_wide_geno_tabular <- function(table, first_snp, last_snp, gds_outfile, rs_pos = 1,
                                   chr_pos = 2, pos_pos = 3){


  valid_IUPAC <- c('A', 'C', 'G', 'T', 'U', 'W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N')
  double_code <- c("AA","TT","CC","GG","AT","TA","AC","CA","AG","GA","TC","CT","TG","GT","CG","GC","NN")


  # read marker columns and transpose to obtain samples x snps
  mt <- t(table[,c(first_snp:last_snp)])
  colnames(mt) <- table[, rs_pos]

  first_row   <- mt[1, -c(1:11)]

  # convert single letter to double
  if (all(first_row %in% valid_IUPAC)) {
    for (call_idx in seq(1:length(valid_IUPAC))) {
      mt[table == valid_IUPAC[call_idx]] <- double_code[call_idx]
    }
  }


  # convert to genind data structure using a separator empty
  gi <- adegenet::df2genind(mt, sep = "", ploidy = 2)
  locna <- gi@loc.n.all
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

  alellelic_dosage <- gi@tab[, ccc]
  alleles <- as.vector(unlist(lapply(gi@all.names,
                                     function(sublist) paste(sublist,
                                                             collapse = "/"))))

  # create gds of snprelate

  SNPRelate::snpgdsCreateGeno(gds.fn = gds_outfile,
                              genmat = t(alellelic_dosage),
                              sample.id = adegenet::indNames(gi),
                              snp.id = adegenet::locNames(gi),
                              snp.chromosome = as.factor(table[,chr_pos]),
                              snp.position = as.integer(table[,pos_pos]),
                              snp.allele = alleles,
                              snpfirstdim = TRUE)
}


read_wide_geno_tabular <- function(table, first_snp, last_snp, gds_outfile, rs_pos = 1,
                                   chr_pos = 2, pos_pos = 3){


  valid_IUPAC <- c('A', 'C', 'G', 'T', 'U', 'W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N')
  double_code <- c("AA","TT","CC","GG","AT","TA","AC","CA","AG","GA","TC","CT","TG","GT","CG","GC","NN")


  # read marker columns and transpose to obtain samples x snps
  mt <- t(table[,c(first_snp:last_snp)])
  colnames(mt) <- table[, rs_pos]

  first_row   <- mt[1, -c(1:11)]

  # convert single letter to double
  if (all(first_row %in% valid_IUPAC)) {
    for (call_idx in seq(1:length(valid_IUPAC))) {
      mt[table == valid_IUPAC[call_idx]] <- double_code[call_idx]
    }
  }


  # convert to genind data structure using a separator empty
  gi <- adegenet::df2genind(mt, sep = "", ploidy = 2)
  locna <- gi@loc.n.all
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

  alellelic_dosage <- gi@tab[, ccc]
  alleles <- as.vector(unlist(lapply(gi@all.names,
                                     function(sublist) paste(sublist,
                                                             collapse = "/"))))

  # create gds of snprelate

  SNPRelate::snpgdsCreateGeno(gds.fn = gds_outfile,
                              genmat = t(alellelic_dosage),
                              sample.id = adegenet::indNames(gi),
                              snp.id = adegenet::locNames(gi),
                              snp.chromosome = as.factor(table[,chr_pos]),
                              snp.position = as.integer(table[,pos_pos]),
                              snp.allele = alleles,
                              snpfirstdim = TRUE)
}



read_longer_geno_tabular <- function(table, gds_outfile) {


}


read_vcf <- function(vcf_path, gds_outfile){
  if (length(vcf_path) > 0){
    SNPRelate::snpgdsVCF2GDS(vcf_path,
                             gds_outfile,
                             method="biallelic.only")
  }
}

read_plink <- function(ped_path, map_path, gds_outfile){

  if ((length(ped_path) > 0) & (length(map_path) > 0)){
    SNPRelate::snpgdsPED2GDS(ped.fn = ped_path,
                             map.fn = map_path,
                             gds_outfile)
  }
}

read_DartR <- function(dart_path,
                       chr_name,
                       pos_name,
                       gds_outfile){

  outfile <- basename(gds_outfile)
  outpath <- dirname(gds_outfile)
  if (length(dart_path) > 0){
    gl <- dartR::gl.read.dart(dart_path)

    # Keep SNPs mapped at nuclear regions
    chr_snps <- as.data.frame(gl@other$loc.metrics) %>%
      mutate(chrom = as.character(!!sym(chr_name)),
             snp_id = gl@loc.names) %>%
      filter(complete.cases(chrom)) %>%
      filter(!!sym(pos_name) > 0)

    chrom_data <- dartR::gl.keep.loc(gl, loc.list = chr_snps$snp_id)

    dartR::gl2gds(chrom_data,
                  snp_pos = pos_name,
                  snp_chr = chr_name,
                  outfile = out_file_gds,
                  outpath = outpath)
  }
}


