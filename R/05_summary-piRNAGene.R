#' Summary method for piRNAGene objects
#'
#' Provides a summary of the piRNA gene information.
#'
#' @param object A piRNAGene object.
#' @return A summary of the piRNA gene information.
#' @name summary-piRNAGene
#' @title Summary Method for piRNAGene Class
#' @aliases summary-piRNAGene summary,piRNAGene-method
#' @rdname summary-piRNAGene
#' @export
#' @importFrom methods callNextMethod
#' @examples
#' pirna_gene <- createPiRNAGene(
#'   8L, "SYMBOL_PI", "piRNA Name",
#'   "piRNA Description", "chr1", 1, 1000,
#'   "+", list(), "pirna1", "PIRNA_SEQ"
#' )
#' summary(pirna_gene)
setMethod("summary", "piRNAGene", function(object) {
  callNextMethod()
  cat("piRNA ID:", getPiRNAID(object), "\n")
  cat("piRNA Sequence:", getPiRNASequence(object), "\n")
})
