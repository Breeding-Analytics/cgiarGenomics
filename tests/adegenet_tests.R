library(adegenet)


gi <- geno$geno
meta <- geno$meta

MAF <- minorAllele(gi)
hist(MAF)

mono <- !isPoly(gi)
length(which(mono))
