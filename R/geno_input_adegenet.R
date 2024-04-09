
read_hapmap <- function(path, first_snp, last_snp, ploidity = 2, sep = ""){
  
  table <- read_tabular_geno(path)
  # read marker columns and transpose to obtain samples x snps
  mt <- t(table[,c(first_snp:last_snp)])
  
  # rs field always in first column for hapmap fmt
  colnames(mt) <- table[, 1]
  
  # convert to genind data structure using a separator empty
  gi <- adegenet::df2genind(mt, sep = sep, ploidy = ploidity)
  
  
  # Format metadata
  
  meta <- data.frame(
    ID = table[,1],
    CHROM = table[,3],
    POS = as.numeric(table[,4])
  )
  out <- list(geno = gi, meta = meta)
  return(out)
}

# Test
# first_snp = 12
# last_snp = 791
# 
# ploidity = 2
# sep = ""
# path <- 'tests/hmp_fmt/geno_pruned_hmp.txt'
# 
# geno <- read_hapmap(path = path,
#                     first_snp = first_snp,
#                     last_snp = last_snp)

read_vcf <- function(path, ploidity){
  vcf <- vcfR::read.vcfR(path)
  mt <- t(vcfR::extract.gt(vcf))
  gi <- adegenet::df2genind(mt, sep = '/', ploidy = ploidity)
  
  meta_vcf <- vcfR::getFIX(vcf)
  
  meta <- data.frame(
    ID = meta_vcf[,'ID'],
    CHROM = meta_vcf[,'CHROM'],
    POS = as.numeric(meta_vcf[,'POS'])
  )
  
  
  out <- list(geno = gi, meta = meta)
  return(out)
}

# Test
#vcf_path <- '~/Projects/2024/Hackathon_BA_Nairobi/cgiarGenomics/tests/vcf_fmt/geno_pruned.vcf.gz'
#geno <- read_vcf(vcf_path, ploidity = 2)

read_DArTSeq_SNP <- function(dart_path, snp_id, chr_name, pos_name){
  
  gl <- dartR::gl.read.dart(dart_path)
  gi <- dartR::gl2gi(gl)
  
  # Prepare metadata
  meta <- data.frame(
    ID = gl@other$loc.metrics[,snp_id],
    CHROM = gl@other$loc.metrics[,chr_name],
    POS = gl@other$loc.metrics[,pos_name]
  )
  
  out <- list(geno = gi, meta = meta)
  return(out)
}

#  Test
# dart_path <- '~/Projects/2024/Hackathon_BA_Nairobi/cgiarGenomics/tests/DartSeq_fmt/DartSeq/Report_DCob23-8528_SNP_mapping_2.csv'
# 
# geno <- read_DArTSeq_SNP(dart_path,
#                          snp_id = "AlleleID",
#                          chr_name = "Chrom_Common_bean_Phaseolus_acutifolius_v10",
#                          pos_name = "ChromPosSnp_Common_bean_Phaseolus_acutifolius_v10")
# 



read_DArTSeq_PA <- function(dart_path, marker_id, chr_name, pos_name){
  
  gl <- dartR::gl.read.silicodart(dart_path)
  gi <- dartR::gl2gi(gl)
  
  # Prepare metadata
  meta <- data.frame(
    ID = gl@other$loc.metrics[,marker_id],
    CHROM = gl@other$loc.metrics[,chr_name],
    POS = gl@other$loc.metrics[,pos_name]
  )
  
  out <- list(geno = gi, meta = meta)
  return(out)
}

# Test
# dart_PA_path <- '~/Projects/2024/Hackathon_BA_Nairobi/cgiarGenomics/tests/DartSeq_fmt/DartSeq/Report_DCob23-8528_SilicoDArT_1.csv'
# 
# geno <- read_DArTSeq_PA(dart_PA_path,
#                         marker_id = "CloneID",
#                         chr_name = "Chrom_Common_bean_Phaseolus_acutifolius_v10",
#                         pos_name = "ChromPosTag_Common_bean_Phaseolus_acutifolius_v10")


read_DArT_Tag <- function(counts.file, dosage.file, ploidy = 4){
  
  # Convert to VCF file
  polyBreedR::dart2vcf(counts.file=counts.file,
                       dosage.file=dosage.file,
                       ploidy = ploidity,
                       vcf.file= "tmp.vcf.gz")
  out <- read_vcf("tmp.vcf.gz", ploidity = ploidity)
  return(out)
}

# test
# counts.file = '~/Projects/2024/Hackathon_BA_Nairobi/cgiarGenomics/tests/DartSeq_fmt/DArT_TAG_sweet_potato/DP23-8340_Allele_match_counts_collapsed.csv'
# dosage.file = '~/Projects/2024/Hackathon_BA_Nairobi/cgiarGenomics/tests/DartSeq_fmt/DArT_TAG_sweet_potato/DP23-8340_Allele_Dose_Report.csv'
#
# ploidity <- 4
# geno <- read_DArT_Tag(counts.file=counts.file,
#                      dosage.file=dosage.file,
#                      ploidy = ploidity)


