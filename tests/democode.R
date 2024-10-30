library(cgiarGenomics)

# Hapmap comparision ------------------------------------------------------

hapMapChar2NumericDouble <- function(hapMap) {
  
  hapMap <- as.data.frame(hapMap)
  dim(hapMap)
  # extract SNP infomation , which is the first 11 columns
  SNPInfo <- hapMap[,1:11]
  
  # remove the first 11 columns
  hapMap <- hapMap[,-c(1:11)]
  missingData=c("NN","FAIL","FAILED","Uncallable","Unused","NA","",-9)
  for(iMiss in missingData){hapMap[which(hapMap==iMiss, arr.ind = TRUE)] <- NA}
  # convert the hapMap to numeric
  Mprov <- t(hapMap); colnames(Mprov) <- SNPInfo[,1]
  hapMapNumeric <- sommer::atcg1234(Mprov, maf = -1, imp = FALSE)
  
  multiAllelic <- setdiff(SNPInfo$`rs#`,colnames(hapMapNumeric$M))
  if(length(multiAllelic) > 0){
    addMulti <- matrix(NA,nrow=nrow(hapMapNumeric$M),ncol=length(multiAllelic))
    addMultiRef <- matrix(NA,nrow=2,ncol=length(multiAllelic))
    colnames(addMulti) <- colnames(addMultiRef) <- multiAllelic
    hapMapNumeric$M <- cbind(hapMapNumeric$M, addMulti)
    hapMapNumeric$ref.alleles <- cbind(hapMapNumeric$ref.alleles, addMultiRef)
  }
  # convert to data frame
  refAlleles <- hapMapNumeric$ref.alleles
  
  hapMapNumeric <- as.data.frame(t(hapMapNumeric$M+1))
  
  # add reference and alternate allele
  SNPInfo$alleles <- apply(refAlleles,2,function(x){paste(na.omit(x),collapse = "/")})
  # convert -9 values to NA
  # hapMapNumeric[hapMapNumeric == -9] <- NA
  
  # get back the column names (accessions)
  colnames(hapMapNumeric) <- colnames(hapMap)
  
  result <- cbind(SNPInfo, hapMapNumeric)
  return(result)
}

hapmap_file <- "~/Projects/2024/Hackathon_BA_Nairobi/cgiarGenomics/tests/hmp_fmt/geno_pruned_hmp.txt"
df <- read_tabular_geno(hapmap_file)

system.time(
  hmp_df <- hapMapChar2NumericDouble(df)
)

system.time(
  gl <- read_hapmap(hapmap_file)
)

# filtering example
gl_maf <- filter_MAF(gl, 0.05)


# Distance 

x <- seploc(gl_maf, n.block=10, parallel=FALSE)
lD <- lapply(x, function(e) dist(as.matrix(e)))
D <- Reduce("+", lD)

library(ape)
library(treeio)
library(ggtree)
library(ggplot2)


tree <- njs(D)

p <- ggtree(tree, layout = 'circular', branch.length='none') +
  geom_tiplab(align=F,
              size=1,
              linetype = 'solid')

system.time(
pca <- adegenet::glPca(gl_maf, nf = 2, useC = T)
)
scatter(pca)


