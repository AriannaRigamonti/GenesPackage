#' Get the tRNA ID
#'
#' This function retrieves the ID of a tRNA gene.
#'
#' @param object A tRNAGene object.
#' @return The tRNA ID.
#' @export
#' @aliases getTRNAID getTRNAID,tRNAGene-method
#' @examples
#' trna_gene <- createTRNAGene(
#'   6L, "SYMBOL_T", "tRNA Name",
#'   "tRNA Description", "chr1", 1, 1000, "+",
#'   list(), "trna1", "TRNA_SEQ"
#' )
#' getTRNAID(trna_gene)
setGeneric("getTRNAID", function(object) standardGeneric("getTRNAID"))

#' @rdname getTRNAID
#' @export
setMethod("getTRNAID", "tRNAGene", function(object) object@tRNAID)
