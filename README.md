# GenesPackage

Author: Arianna Rigamonti

## Description
**GenesPackage** is an R package that provides a structured representation of different types of genes. It includes methods for creating, manipulating, and accessing information related to various gene classes, such as protein-coding genes, lncRNA genes, and microRNA genes.

This package is developed as part of the **Scientific Programming** course project.

## Features
- **Creation of gene objects**: Define different gene types, including lncRNAs, microRNAs, piRNAs, and more.
- **Gene information retrieval**: Access attributes like gene symbol, description, chromosome location, and product.
- **Gene modification**: Update attributes dynamically within an R session.
- **Integration with GenomicRanges**: Utilize `GRanges` objects for genomic coordinate representation.

## Installation
To install the package from source, use the following command:

```r
# Install from local tar.gz file
install.packages("GenesPackage_0.1.0.tar.gz", repos = NULL, type = "source")

# Load the package
library(GenesPackage)
```

## Usage

1. Creating a Generic Gene Object

```r
library(GenomicRanges)

# Define genomic coordinates
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1, end = 1000))

# Create a generic gene object
gene <- createGene(ID = 1L, symbol = "SYMBOL", name = "Gene Name",
                   description = "Description", chr = "chr1", start = 1, 
                   end = 1000, strand = "+", product = list())

# Display gene information
gene
```
2. Creating an lncRNA Gene Object

```r
lncrna_gene <- createLncRNAGene(ID = 2L, symbol = "SYMBOL_LNC", name = "LncRNA Name",
                                description = "LncRNA Description", chr = "chr1", 
                                start = 1, end = 1000, strand = "+", 
                                product = list(), lncRNAID = "lncrna1", 
                                RNASequence = "RNA_SEQUENCE")

lncrna_gene
```

3. Retrieving Gene Attributes
   
```r
getID(gene)        # Get gene ID
getSymbol(gene)    # Get gene symbol
getName(gene)      # Get gene name
getDescription(gene) # Get gene description
getStructure(gene) # Get genomic coordinates
```

4. Modifying Gene Attributes
   
```r  
setSymbol(gene) <- "NEW_SYMBOL"
setName(gene) <- "New Gene Name"
setDescription(gene) <- "Updated Description"
```

## Dependencies
	•	R version: 4.4.1 or later
	•	Required packages: GenomicRanges, IRanges, S4Vectors, BiocGenerics

To install dependencies:

```r
install.packages(c("GenomicRanges", "IRanges", "S4Vectors", "BiocGenerics"), dependencies = TRUE)
```
