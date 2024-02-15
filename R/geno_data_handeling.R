
# Read functions ----------------------------------------------------------

read_wide_geno_tabular <- function(table, first_snp, last_snp, gds_outfile, rs_pos = 1,
                                   chr_pos = 2, pos_pos = 3){

  # read marker columns and transpose to obatin samples x snps
  mt <- t(table[,c(first_snp:last_snp)])
  colnames(mt) <- table[, rs_pos]
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

read_DartR <- function(dart_path, gds_outfile){
  outfile <- basename(gds_outfile)
  outpath <- dirname(gds_outfile)
  if (length(dart_path) > 0){
    genofile <- dartR::gl.read.dart(dart_path)
    dartR::gl2gds(genofile,
           outfile = outfile,
           outpath = outpath)
  }
}


