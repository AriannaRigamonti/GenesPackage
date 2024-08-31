#' Summary method for lncRNAGene objects
#'
#' Provides a summary of the lncRNA gene information.
#'
#' @param object A lncRNAGene object.
#' @return A summary of the lncRNA gene information.
#' @name summary-lncRNAGene
#' @title Summary Method for lncRNAGene Class
#' @aliases summary-lncRNAGene summary,lncRNAGene-method
#' @rdname summary-lncRNAGene
#' @export
#' @importFrom methods callNextMethod
#' @examples
#' lncrna_gene <- createLncRNAGene(
#'   2L, "SYMBOL_LNC", "LncRNA Name",
#'   "LncRNA Description", "chr1", 1, 1000,
#'   "+", list(), "lncrna1", "RNA_SEQUENCE"
#' )
#' summary(lncrna_gene)
setMethod("summary", "lncRNAGene", function(object) {
  callNextMethod()
  cat("lncRNA ID:", getLncRNAID(object), "\n")
  cat("RNA Sequence:", getRNASequence(object), "\n")
})
