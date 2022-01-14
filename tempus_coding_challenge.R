### Tempus coding challenge

##Load packages
library(VariantAnnotation)
library(SeqVarTools)
library(here)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

##read vcf file into GRanges object 
vcf <- readVcf(here::here("Challenge_data_(1).vcf"), "hg19")

## 0. set up annotation dataframe
ann <- as.data.frame(ranges(vcf))

## 1. Type of variation (substitution, insertion, CNV, etc.) and their effect (missense, silent, 
# intergenic, etc.). If there are multiple effects, annotate with the most deleterious 
# possibility.

#Determine type of variation with SeqVarTools package
ann$var.type[which(isSubstitution(vcf)==TRUE)] <- "substitution"
ann$var.type[which(isInsertion(vcf)==TRUE)] <- "insertion"  #no insertions???
ann$var.type[which(isDeletion(vcf)==TRUE)] <- "deletion" #no deletions???

#Effect of variation (missence (changes the amino acid), silent (doesn't change amino acid), intergenic)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(vcf)[1:25] <- paste0("chr",seqlevels(vcf))
common.seq.levels <- intersect(seqlevels(vcf), seqlevels(txdb))
vcf <- keepSeqlevels(vcf, common.seq.levels, pruning.mode="coarse")
txbd <- keepSeqlevels(txdb, common.seq.levels, pruning.mode="coarse")

coding <- predictCoding(vcf, txdb, seqSource=Hsapiens) #compute amino acid changes

ann$frameshift <- FALSE
ann$missence <- FALSE
ann$silent <- FALSE
  
for (i in 1:length(ranges(vcf))){
  indices = which(ranges(coding) == ranges(vcf)[i])
  if(length(indices) == 0){
    next
  }
  for (j in 1:length(indices)){
    index = indices[j] #index = indices[1]
    if(toString(coding$VARAA[index])==""){
      ann$frameshift[i] = TRUE
    }
    else if (toString(coding$VARAA[index])!=toString(coding$REFAA[index])){
      ann$missence[i] =TRUE
    }
    else if (toString(coding$VARAA[index])==toString(coding$REFAA[index])){
      ann$silent[i] =TRUE
    }
    
  }
  
}


#https://bioconductor.org/packages/release/bioc/vignettes/VariantAnnotation/inst/doc/VariantAnnotation.pdf
#proteinCoding()

## 2. Depth of sequence coverage at the site of variation
ann$depth <- info(vcf)$DP #pulls out depth from info field DP

## 3. Number of reads supporting the variant.

# DRPA is this ratio of depths such that ALT Depth:REF Depth and is presented as a decimal so
# DPRA = ALT Depth / REF Depth for a given variant. To find the number of reads supporting the ALT variant we must 
# solve the following equation DPRA = ALT / (DP - ALT) for ALT. ALT = (DPRA * DP) / (1 + DPRA)
dpra <- numeric()
for(i in 1:6977){
  dpra[i] <- info(vcf)$DPRA[[i]][1]
}
dpra <- unlist(dpra)
ann$var.num <- dpra*info(vcf)$DP / (1 + dpra)

## 4. Percentage of reads supporting the variant versus those supporting reference reads. 

# percentage of reads supporting Variant
ann$percent.var <- ann$var.num/ann$depth

#percentage of reads supporting ref
ann$percent.ref <- (ann$depth - ann$var.num)/ann$depth

## 5. 

## 6. 


write.csv(ann,here::here(), row.names = TRUE)



