# Tempus Coding Challenge

Script to annotated genetic variants from VCF file.

## Description

This script first takes a .vcf file, converts it a GRanges data object. It then creates a separate data frame of annotations including the depth of sequence coverage, variant type, alternate allele frequency, and a few others. Lastly, it saves this annotation data frame as a .csv file in the working directory. 

## Getting Started

### Dependencies

* Make sure you have the following packages installed from CRAN or Bioconductor.

```
BiocManager::install("VariantAnnotation")
BiocManager::install("SeqVarTools")
install.packages("here")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGen")
```

### Installing

* No installation necessary

### Executing program

* Create a directory that contains the script and the VCF data file. 
* Open R from that directory.
* Simply run script from top to bottom. 

## Help

Contact the Author.

## Authors

ex. Abby Spangler  
ex. [LinkedIn](www.linkedin.com/in/abby-spangler-72166883)

## Acknowledgments

Inspiration, most helpful packages, etc.
* [VariantAnnotation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html)
* [SeqVarTools](http://www.bioconductor.org/packages/release/bioc/html/SeqVarTools.html)

Thank you for interviewing me, Tempus!
