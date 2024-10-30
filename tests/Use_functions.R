library(dplyr)
library(stringr)

# READ TABULAR hapmap like ------------------------------------------------

hmp_path <- "~/Projects/2024/Hackathon_BA_Nairobi/cgiarGenomics/tests/hmp_fmt/geno_pruned_hmp.txt"
read_hapmap(hmp_path)


# VCF reading -------------------------------------------------------------

vcf_path <-"~/Projects/2024/Hackathon_BA_Nairobi/cgiarGenomics/tests/vcf_fmt/geno_pruned.vcf.gz"
read_vcf(vcf_path)


# DartSeq Reading ---------------------------------------------------------

dart_seq_path <- "~/Projects/2024/Hackathon_BA_Nairobi/cgiarGenomics/tests/DartSeq_fmt/DartSeq/Report_DCob23-8528_SNP_mapping_2.csv"
chr_col <- "Chrom_Common_bean_Phaseolus_acutifolius_v10"
pos_col <- "ChromPosTag_Common_bean_Phaseolus_acutifolius_v10"
snp_id_col <- "AlleleID"
read_DArTSeq_SNP(dart_seq_path,
                 snp_id_col,
                 chr_col,
                 pos_col)


# DartTag -----------------------------------------------------------------

dartTag_2n_path_snp <- "~/Projects/2024/Hackathon_BA_Nairobi/cgiarGenomics/tests/DartSeq_fmt/DArT_TAG_bean/Report_DCob19-4728_tgSNP_1.csv"
dartTag_2n_path_count <- "~/Projects/2024/Hackathon_BA_Nairobi/cgiarGenomics/tests/DartSeq_fmt/DArT_TAG_bean/Report_DCob19-4728_tgRawCounts_2.csv"

read_DArT_Tag(dartTag_2n_path_count,
              dartTag_2n_path_snp,
              ploidity = 2)

polyBreedR::dart2vcf(dartTag_2n_path_count, dartTag_2n_path_snp, 'test.vcf', 2, 7)


# Test dartag reading function --------------------------------------------

dartTag_2n_path_snp <- "~/Projects/2024/Hackathon_BA_Nairobi/cgiarGenomics/tests/DartSeq_fmt/DArT_TAG_sweet_potato/DP23-8340_Allele_Dose_Report.csv" 
dartTag_2n_path_counts <- "~/Projects/2024/Hackathon_BA_Nairobi/cgiarGenomics/tests/DartSeq_fmt/DArT_TAG_sweet_potato/DP23-8340_Allele_match_counts_collapsed.csv" 
data <- read.csv(dartTag_2n_path_snp, header=F, check.names=F)
cols <- 6:ncol(data)
first.data.row <- 9

# rows with data
rows <- first.data.row:nrow(data)
# sample names
id <- as.character(data[first.data.row-2,cols])
# Samples with duplicated id?
dupes <- unique(id[which(duplicated(id))])

n <- length(id)
# sub df with snp_id, chr, pos
map <- data[rows,c(1,4,5)]
m <- nrow(map)
colnames(map) <- c("marker","chrom","position")
map$position <- as.integer(map$position)
ploidy <- 4
geno <- ploidy - apply(data[rows,cols],2,as.integer)
dimnames(geno) <- list(map$marker,id)

GT <- apply(geno, c(1, 2), make_GT, ploidy = ploidy)
data2 <- read.csv(counts.file, header = F, check.names = F)
cols <- 6:ncol(data2)
rows <- first.data.row:nrow(data2)
id <- as.character(data2[first.data.row - 2, cols])


data2 <- read.csv(dartTag_2n_path_counts, header = F, check.names = F)
cols <- 6:ncol(data2)
rows <- first.data.row:nrow(data2)
id <- as.character(data2[first.data.row - 2, cols])
if (length(dupes) > 0) {
  ix <- which(id %in% dupes)
  id[ix] <- apply(cbind(id[ix], as.character(data2[first.data.row - 
                                                     1, cols])[ix]), 1, paste, collapse = ".")
}
allele.seq <- matrix(data2[rows, 3], ncol = 2, byrow = T)

REF.ALT <- t(apply(allele.seq, 1, function(x) {
  bases <- c("A", "C", "G", "T")
  hap.ref <- unlist(strsplit(x[1], split = ""))
  hap.alt <- unlist(strsplit(x[2], split = ""))
  k <- which(hap.ref != hap.alt & hap.ref %in% bases & 
               hap.alt %in% bases)
  if (length(k) == 1) {
    return(c(hap.ref[k], hap.alt[k]))
  }
  else {
    return(c("N", "N"))
  }
}))

make_GT <- function (x, ploidy) 
{
  if (is.na(x)) {
    y <- rep(".", ploidy)
  }
  else {
    y <- c(rep(0, ploidy - x), rep(1, x))
  }
  return(paste(y, collapse = "/"))
}


# read dart ---------------------------------------------------------------

get_dart_gt_data(dartTag_2n_path_snp,
                sample.name.row = 7,
                first.data.row = 9,
                first.data.col = 6,
                snp_id_col = 'MarkerID')

dartTag_2n_path_snp <- "~/Projects/2024/Hackathon_BA_Nairobi/cgiarGenomics/tests/DartSeq_fmt/DArT_TAG_bean/DCob20-5507_Standard_DArT_report.csv" 
x <- read.csv(dartTag_2n_path_snp, header = F)
meta_data_test <- x[8:nrow(x),1:4]
colnames(meta_data_test) <- x[7,1:4]

purrr::map2_dfr(x[8:nrow(x),2],
                x[8:nrow(x),3],\(.x,.y) compare_tags(.x, .y))

datr_2row_path <- "~/Projects/2024/Hackathon_BA_Nairobi/cgiarGenomics/tests/DartSeq_fmt/DArT_TAG_bean/Report_DCob20-5507_2_row.csv" 
x <- read.csv(datr_2row_path, header = F)
meta_data_test <- x[8:nrow(x),1:4]
colnames(meta_data_test) <- x[7,1:4]

tidyr::pivot_wider(meta_data_test, )
  
pivot_
meta_data_test
