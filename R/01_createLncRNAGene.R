#' Create an lncRNAGene object
#'
#' @param ID integer The gene ID (e.g., Ensembl ID or NCBI gene ID).
#' @param symbol character. The gene symbol.
#' @param name character. The gene name.
#' @param description character. The gene description.
#' @param chr character. The chromosome location.
#' @param start numeric. The start position of the gene.
#' @param end numeric. The end position of the gene.
#' @param strand character. The strand of the gene.
#' @param product list. The gene product.
#' @param lncRNAID character. The lncRNA ID.
#' @param RNASequence character. The RNA sequence.
#' @return An lncRNAGene object.
#' @export
#' @examples
#' lncrna_gene <- createLncRNAGene(
#'    2L, "SYMBOL_LNC", "LncRNA Name",
#'   "LncRNA Description", "chr1", 1, 1000, "+",
#'   list(), "lncrna1", "RNA_SEQUENCE"
#' )
#' lncrna_gene
createLncRNAGene <- function(ID, symbol, name, description, chr, start, end, strand,
                             product, lncRNAID, RNASequence) {
  structure <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = start, end = end),
    strand = strand
  )
  new("lncRNAGene",
    ID = ID, symbol = symbol, name = name, description = description,
    structure = structure, product = product, lncRNAID = lncRNAID, RNASequence = RNASequence
  )
}
