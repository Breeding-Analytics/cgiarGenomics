library(tidyverse)
library(gdsfmt)
library(ggrepel)
library(ggpubr)
library(vcfR)
setwd("~/Projects/2024/Hackathon_BA_Nairobi/cgiarGenomics")

hmp_path <- 'tests/hmp_fmt/geno_pruned_hmp.txt'
gds_outfile <- 'tests/hmp_fmt/geno_hmp.gds'


first_snp = 12
last_snp = 791
chr_pos = 2
pos_pos = 3
ploidity = 2
rs_pos = 1
sep = ""

missingData=c("NN","FAIL","FAILED","Uncallable","Unused","NA","",-9)

df <- as.data.frame(data.table::fread(hmp_path,
                                      sep = '\t',
                                      header = TRUE,
                                      na.strings = missingData))
out_file_gds <- "./tests/hmp_geno_tetraploid.gds"


table <- df


# read marker columns and transpose to obtain samples x snps
mt <- t(table[,c(first_snp:last_snp)])
colnames(mt) <- table[, rs_pos]

# convert to genind data structure using a separator empty
gi <- adegenet::df2genind(mt, sep = sep, ploidy = ploidity)
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
alellelic_dosage <- alellelic_dosage +
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


genofile <- SNPRelate::snpgdsOpen(gds_outfile)

plot_missingnes_by_marker(genofile)
plot_overall_missingness(genofile)
plot_MAF(genofile)


SNPRelate::snpgdsSNPRateFreq(genofile)

gt <- SNPRelate::snpgdsGetGeno(genofile)

read_vcf('tests/DartSeq_fmt/DArTag_poly.vcf.gz', )

SNPRelate::snpgdsSNPRateFreq(genofile)

genofile <- SNPRelate::snpgdsOpen('tests/DartSeq_fmt/DArTag_poly.gds')
plot_missingnes_by_marker(genofile)


vcf <- read.vcfR('tests/DartSeq_fmt/DArTag_poly.vcf.gz')

vcf@gt[1:5,1:5]
