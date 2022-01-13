#tempus coding challenge

library(VariantAnnotation)
library(SeqVarTools)
library(here)

library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#read vcf file into GRanges object 
vcf <- readVcf(here::here("Challenge_data_(1).vcf"), "hg19")

info(vcf)

# 0. set up annotation dataframe
ann <- as.data.frame(rowRanges(vcf))

# 1. Type of variation (substitution, insertion, CNV, etc.) and their effect (missense, silent, 
# intergenic, etc.). If there are multiple effects, annotate with the most deleterious 
# possibility.
ann$var.type[which(isSubstitution(vcf)==TRUE)] <- "substitution"
ann$var.type[which(isInsertion(vcf)==TRUE)] <- "insertion"  #no insertions???
ann$var.type[which(isDeletion(vcf)==TRUE)] <- "deletion" #no deletions???

#https://bioconductor.org/packages/release/bioc/vignettes/VariantAnnotation/inst/doc/VariantAnnotation.pdf
#proteinCoding()

# 2. Depth of sequence coverage at the site of variation
ann$depth <- info(vcf)$DP #pulls out depth from info field DP

# 3. Number of reads supporting the variant.
# My understaning is taht DRPA is this ratio of depths such that ALT Depth:REF Depth and is presented as a decimal so
# DPRA = ALT Depth / REF Depth for a given variant. To find the number of reads supporting the ALT variant we must 
# solve the following equation DPRA = ALT / (DP - ALT) for ALT. ALT = (DPRA * DP) / (1 + DPRA)

df.dpra <- numeric()
for(i in 1:6977){
  df.dpra[i] <- info(vcf)$DPRA[[i]][1]
}
df.dpra <- unlist(df.dpra)
ann$var.num <- df.dpra*info(vcf)$DP / (1 + df.dpra) #DPRA from info field. multiply this number by depth to see how many reads called the alternative 

# 4. Percentage of reads supporting the variant versus those supporting reference reads. 

# percentage of reads supporting Variant
ann$percent.var <- ann$var.num/ann$depth

#percentage of reads supporting ref
ann$percent.red <- (ann$depth - ann$var.num)/ann$depth


colData(vcf)
rowData(vcf)
rowRanges(vcf) #similar to rownames



