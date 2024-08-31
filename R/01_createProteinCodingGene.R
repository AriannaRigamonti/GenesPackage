#' Create a ProteinCodingGene object
#'
#' This function creates a new ProteinCodingGene object.
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
#' @param proteinID character. The protein ID.
#' @param proteinSequence character. The protein sequence.
#' @return A ProteinCodingGene object.
#' @export
#' @examples
#' protein_gene <- createProteinCodingGene(
#'   1L, "SYMBOL", "Gene Name", "Description",
#'   "chr1", 1, 1000, "+", list(), "protein1", "SEQUENCE"
#' )
#' protein_gene
createProteinCodingGene <- function(ID, symbol, name, description, chr, start, end, strand,
                                    product, proteinID, proteinSequence) {
  structure <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = start, end = end),
    strand = strand
  )
  new("ProteinCodingGene",
    ID = ID, symbol = symbol, name = name, description = description,
    structure = structure, product = product, proteinID = proteinID,
    proteinSequence = proteinSequence
  )
}
