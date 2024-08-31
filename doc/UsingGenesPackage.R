## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
# Install the package from source
# install.packages("GenesPackage_0.1.0.tar.gz", repos = NULL, type = "source")

# Load the package
library(GenesPackage)

## -----------------------------------------------------------------------------
library(GenomicRanges)
# Create a GRanges object for the gene structure
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1, end = 1000))

# Create a generic Gene object
gene <- createGene(ID = 1L, symbol = "SYMBOL", name = "Gene Name", 
                   description = "Description", chr = "chr1", start = 1, 
                   end = 1000, strand = "+", product = list())

# Display the gene object
gene

## -----------------------------------------------------------------------------
# Create a lncRNAGene object
lncrna_gene <- createLncRNAGene(ID = 2L, symbol = "SYMBOL_LNC", 
                                name = "LncRNA Name", description = "LncRNA Description", 
                                chr = "chr1", start = 1, end = 1000, strand = "+", 
                                product = list(), lncRNAID = "lncrna1", 
                                RNASequence = "RNA_SEQUENCE")

# Display the lncRNA gene object
lncrna_gene

## -----------------------------------------------------------------------------
# Create a microRNAGene object
mirna_gene <- createMicroRNAGene(ID = 3L, symbol = "SYMBOL_MIR", 
                                 name = "MicroRNA Name", description = "MicroRNA Description", 
                                 chr = "chr1", start = 1, end = 1000, strand = "+", 
                                 product = list(), microRNAID = "mirna1", 
                                 microRNASequence = "SEED_SEQ")

# Display the microRNA gene object
mirna_gene

## -----------------------------------------------------------------------------
# Create a piRNAGene object
pirna_gene <- createPiRNAGene(ID = 8L, symbol = "SYMBOL_PI", 
                              name = "piRNA Name", description = "piRNA Description", 
                              chr = "chr1", start = 1, end = 1000, strand = "+", 
                              product = list(), piRNAID = "pirna1", 
                              piRNASequence = "PIRNA_SEQ")

# Display the piRNA gene object
pirna_gene

## -----------------------------------------------------------------------------
# Create a ProteinCodingGene object
protein_gene <- createProteinCodingGene(ID = 1L, symbol = "SYMBOL", 
                                        name = "Gene Name", description = "Description", 
                                        chr = "chr1", start = 1, end = 1000, strand = "+", 
                                        product = list(), proteinID = "protein1", 
                                        proteinSequence = "SEQUENCE")

# Display the protein-coding gene object
protein_gene

## -----------------------------------------------------------------------------
# Create an rRNAGene object
rrna_gene <- createRRNAGene(ID = 7L, symbol = "SYMBOL_R", 
                            name = "rRNA Name", description = "rRNA Description", 
                            chr = "chr1", start = 1, end = 1000, strand = "+", 
                            product = list(), rRNAID = "rrna1", rRNASequence = "RRNA_SEQ")

# Display the rRNA gene object
rrna_gene

## -----------------------------------------------------------------------------
# Create a siRNAGene object
sirna_gene <- createSiRNAGene(ID = 4L, symbol = "SYMBOL_SI", 
                              name = "siRNA Name", description = "siRNA Description", 
                              chr = "chr1", start = 1, end = 1000, strand = "+", 
                              product = list(), siRNAID = "sirna1", 
                              siRNASequence = "SIRNA_SEQ")

# Display the siRNA gene object
sirna_gene

## -----------------------------------------------------------------------------
# Create a snoRNAGene object
snorna_gene <- createSnoRNAGene(ID = 5L, symbol = "SYMBOL_SNO", 
                                name = "snoRNA Name", description = "snoRNA Description", 
                                chr = "chr1", start = 1, end = 1000, strand = "+", 
                                product = list(), snoRNAID = "snorna1", 
                                snoRNASequence = "SNORNA_SEQ")

# Display the snoRNA gene object
snorna_gene

## -----------------------------------------------------------------------------
# Create a tRNAGene object
trna_gene <- createTRNAGene(ID = 6L, symbol = "SYMBOL_T", 
                            name = "tRNA Name", description = "tRNA Description", 
                            chr = "chr1", start = 1, end = 1000, strand = "+", 
                            product = list(), tRNAID = "trna1", 
                            tRNASequence = "TRNA_SEQ")

# Display the tRNA gene object
trna_gene

## -----------------------------------------------------------------------------
# Get the gene ID
getID(gene)

## -----------------------------------------------------------------------------
# Get the gene symbol
getSymbol(gene)

## -----------------------------------------------------------------------------
# Get the gene name
getName(gene)

## -----------------------------------------------------------------------------
# Get the gene description
getDescription(gene)

## -----------------------------------------------------------------------------
# Get the gene structure
getStructure(gene)

## -----------------------------------------------------------------------------
# Get the gene product
getProduct(gene)

## -----------------------------------------------------------------------------
# Set a new gene ID
setID(gene) <- 2L
getID(gene)

## -----------------------------------------------------------------------------
# Set a new gene symbol
setSymbol(gene) <- "NEW_SYMBOL"
getSymbol(gene)

## -----------------------------------------------------------------------------
# Set a new gene name
setName(gene) <- "New Gene Name"
getName(gene)

## -----------------------------------------------------------------------------
# Set a new gene description
setDescription(gene) <- "New Description"
getDescription(gene)

## -----------------------------------------------------------------------------
# Set a new gene structure
new_structure <- GRanges(seqnames = "chr2", ranges = IRanges(start = 100, end = 2000), strand = "-")
setStructure(gene) <- new_structure
getStructure(gene)

## -----------------------------------------------------------------------------
# Set a new gene product
setProduct(gene) <- list("new product")
getProduct(gene)

## -----------------------------------------------------------------------------
# Summary of a generic Gene object
summary(gene)

## -----------------------------------------------------------------------------
# Summary of a lncRNA Gene object
summary(lncrna_gene)

## -----------------------------------------------------------------------------
# Summary of a microRNA Gene object
summary(mirna_gene)

## -----------------------------------------------------------------------------
# Summary of a piRNA Gene object
summary(pirna_gene)

## -----------------------------------------------------------------------------
# Summary of a Protein-Coding Gene object
summary(protein_gene)

## -----------------------------------------------------------------------------
# Summary of an rRNA Gene object
summary(rrna_gene)

## -----------------------------------------------------------------------------
# Summary of a siRNA Gene object
summary(sirna_gene)

## -----------------------------------------------------------------------------
# Summary of a snoRNA Gene object
summary(snorna_gene)

## -----------------------------------------------------------------------------
# Summary of a tRNA Gene object
summary(trna_gene)

## -----------------------------------------------------------------------------
sessionInfo()

