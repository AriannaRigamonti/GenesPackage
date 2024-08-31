#' Get the snoRNA sequence
#'
#' This function retrieves the sequence of a snoRNA gene.
#'
#' @param object A snoRNAGene object.
#' @return The snoRNA sequence.
#' @export
#' @aliases getSnoRNASequence getSnoRNASequence,snoRNAGene-method
#' @examples
#' snorna_gene <- createSnoRNAGene(
#'   5L, "SYMBOL_SNO", "snoRNA Name",
#'   "snoRNA Description", "chr1", 1, 1000, "+",
#'   list(), "snorna1", "SNORNA_SEQ"
#' )
#' getSnoRNASequence(snorna_gene)
setGeneric("getSnoRNASequence", function(object) standardGeneric("getSnoRNASequence"))

#' @rdname getSnoRNASequence
#' @export
setMethod("getSnoRNASequence", "snoRNAGene", function(object) object@snoRNASequence)
