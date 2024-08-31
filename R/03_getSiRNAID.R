#' Get the siRNA ID
#'
#' This function retrieves the ID of a siRNA gene.
#'
#' @param object A siRNAGene object.
#' @return The siRNA ID.
#' @export
#' @aliases getSiRNAID getSiRNAID,siRNAGene-method
#' @examples
#' sirna_gene <- createSiRNAGene(
#'   4L, "SYMBOL_SI", "siRNA Name",
#'   "siRNA Description", "chr1", 1, 1000, "+",
#'   list(), "sirna1", "SIRNA_SEQ"
#' )
#' getSiRNAID(sirna_gene)
setGeneric("getSiRNAID", function(object) standardGeneric("getSiRNAID"))

#' @rdname getSiRNAID
#' @export
setMethod("getSiRNAID", "siRNAGene", function(object) object@siRNAID)
