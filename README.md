# Tempus Coding Challenge

Script to annotated genetic variants from VCF file.

## Description

This script first takes a .vcf file and converts it a GRanges data object. It then creates a separate data frame of annotations including the depth of sequence coverage, variant type, alternate allele frequency, and a few others. Lastly, it saves this annotation data frame as a .csv file in the working directory. 

## Locations of answers to challenge questions in output file
1. `var.type` and `var.effect` columns
2. `depth` column
3. `var.num` column
4. `percent.var` and `percent.ref` columns
5. not completed
6. `perc.ref.for`

[Link to GitHub repo](https://github.com/abspangler13/tempus_coding_challenge)

## Getting Started

### Dependencies

* Make sure you have the following packages installed from CRAN or Bioconductor.

```
BiocManager::install("VariantAnnotation")
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

Thank you for interviewing me, Tempus!
