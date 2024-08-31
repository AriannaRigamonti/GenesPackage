#' Summary method for ProteinCodingGene objects
#'
#' Provides a summary of the protein-coding gene information.
#'
#' @param object A ProteinCodingGene object.
#' @return A summary of the protein-coding gene information.
#' @name summary-ProteinCodingGene
#' @export
#' @aliases summary-ProteinCodingGene summary,ProteinCodingGene-method
#' @importFrom methods callNextMethod
#' @examples
#' protein_gene <- createProteinCodingGene(
#'   1L, "SYMBOL", "Gene Name",
#'   "Description", "chr1", 1, 1000, "+",
#'   list(), "protein1", "SEQUENCE"
#' )
#' summary(protein_gene)
setMethod("summary", "ProteinCodingGene", function(object) {
  callNextMethod()
  cat("Protein ID:", getProteinID(object), "\n")
  cat("Protein Sequence:", getProteinSequence(object), "\n")
})
