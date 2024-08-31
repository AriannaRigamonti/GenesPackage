#' Get the siRNA sequence
#'
#' This function retrieves the sequence of a siRNA gene.
#'
#' @param object A siRNAGene object.
#' @return The siRNA sequence.
#' @export
#' @aliases getSiRNASequence getSiRNASequence,siRNAGene-method
#' @examples
#' sirna_gene <- createSiRNAGene(
#'   4L, "SYMBOL_SI", "siRNA Name",
#'   "siRNA Description", "chr1", 1, 1000, "+",
#'   list(), "sirna1", "SIRNA_SEQ"
#' )
#' getSiRNASequence(sirna_gene)
setGeneric("getSiRNASequence", function(object) standardGeneric("getSiRNASequence"))

#' @rdname getSiRNASequence
#' @export
setMethod("getSiRNASequence", "siRNAGene", function(object) object@siRNASequence)
