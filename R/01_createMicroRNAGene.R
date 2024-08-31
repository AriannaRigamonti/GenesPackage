#' Create a microRNAGene object
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
#' @param microRNAID character. The microRNA ID.
#' @param microRNASequence character. The microRNA sequence.
#' @return A microRNAGene object.
#' @export
#' @examples
#' mirna_gene <- createMicroRNAGene(
#'   3L, "SYMBOL_MIR", "MicroRNA Name",
#'   "MicroRNA Description", "chr1", 1, 1000, "+",
#'   list(), "mirna1", "SEED_SEQ"
#' )
#' mirna_gene
createMicroRNAGene <- function(ID, symbol, name, description, chr, start, end, strand, product, microRNAID, microRNASequence) {
  structure <- GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = start, end = end), strand = strand)
  new("microRNAGene", ID = ID, symbol = symbol, name = name, description = description, structure = structure, product = product, microRNAID = microRNAID, microRNASequence = microRNASequence)
}
