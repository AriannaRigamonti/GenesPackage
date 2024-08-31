#' Get the piRNA sequence
#'
#' This function retrieves the sequence of a piRNA gene.
#'
#' @param object A piRNAGene object.
#' @return The piRNA sequence.
#' @export
#' @aliases getPiRNASequence getPiRNASequence,piRNAGene-method
#' @examples
#' pirna_gene <- createPiRNAGene(
#'   8L, "SYMBOL_PI", "piRNA Name",
#'   "piRNA Description", "chr1", 1, 1000, "+",
#'   list(), "pirna1", "PIRNA_SEQ"
#' )
#' getPiRNASequence(pirna_gene)
setGeneric("getPiRNASequence", function(object) standardGeneric("getPiRNASequence"))

#' @rdname getPiRNASequence
#' @export
setMethod("getPiRNASequence", "piRNAGene", function(object) object@piRNASequence)
