#' Create a Gene object
#'
#' This function creates a new Gene object.
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
#' @return A Gene object.
#' @export
#' @examples
#' gene <- createGene(1L, "SYMBOL", "Gene Name", "Description", "chr1", 1, 1000, "+", list())
#' gene
createGene <- function(ID, symbol, name, description, chr, start, end, strand, product) {
  structure <- GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = start, end = end), strand = strand)
  new("Gene", ID = ID, symbol = symbol, name = name, description = description, structure = structure, product = product)
}
