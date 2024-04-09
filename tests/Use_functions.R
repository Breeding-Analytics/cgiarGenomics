library(gdsfmt)
library(adegenet)
library(SNPRelate)
library(SeqArray)
library(dplyr)
library(stringr)


# READ TABULAR hapmap like ------------------------------------------------

hmp_path <- "~/Projects/2024/Hackathon_BA_Nairobi/cgiarGenomics/tests/hmp_fmt/geno_pruned_hmp.txt"

missingData=c("NN","FAIL","FAILED","Uncallable","Unused","NA","",-9)

df <- as.data.frame(data.table::fread(hmp_path,
                                      sep = '\t',
                                      header = TRUE,
                                      na.strings = missingData))
out_file_gds <- "./tests/hmp_geno_tetraploid.gds"
gi <- read_wide_geno_tabular(table = df,
                             first_snp = 12,
                             last_snp = 791,
                             rs_pos = 1,
                             chr_pos = 3,
                             pos_pos = 4,
                             gds_outfile = out_file_gds,
                             ploidity = 4,
)


genofile <- SNPRelate::snpgdsOpen(out_file_gds)

per_marker <- SNPRelate::snpgdsSNPRateFreq(genofile)
per_sample <- SNPRelate::snpgdsSampMissRate(genofile)

SNPRelate::snpgdsSNPList(genofile)



read.gdsn(index.gdsn(genofile, "genotype"))



SNPRelate::snpgdsClose(genofile)

