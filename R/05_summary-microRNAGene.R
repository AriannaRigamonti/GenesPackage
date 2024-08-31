#' Summary method for microRNAGene objects
#'
#' Provides a summary of the microRNA gene information.
#'
#' @param object A microRNAGene object.
#' @return A summary of the microRNA gene information.
#' @name summary-microRNAGene
#' @title Summary Method for microRNAGene Class
#' @aliases summary-microRNAGene summary,microRNAGene-method
#' @rdname summary-microRNAGene
#' @export
#' @importFrom methods callNextMethod
#' @examples
#' mirna_gene <- createMicroRNAGene(
#'   3L, "SYMBOL_MIR", "MicroRNA Name",
#'   "MicroRNA Description", "chr1", 1, 1000,
#'   "+", list(), "mirna1", "SEED_SEQ"
#' )
#' summary(mirna_gene)
setMethod("summary", "microRNAGene", function(object) {
  callNextMethod()
  cat("microRNA ID:", getMicroRNAID(object), "\n")
  cat("microRNA Sequence:", getMicroRNASequence(object), "\n")
})
