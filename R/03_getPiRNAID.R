#' Get the piRNA ID
#'
#' This function retrieves the ID of a piRNA gene.
#'
#' @param object A piRNAGene object.
#' @return The piRNA ID.
#' @export
#' @aliases getPiRNAID getPiRNAID,piRNAGene-method
#' @examples
#' pirna_gene <- createPiRNAGene(
#'   8L, "SYMBOL_PI", "piRNA Name",
#'   "piRNA Description", "chr1", 1, 1000, "+",
#'   list(), "pirna1", "PIRNA_SEQ"
#' )
#' getPiRNAID(pirna_gene)
setGeneric("getPiRNAID", function(object) standardGeneric("getPiRNAID"))

#' @rdname getPiRNAID
#' @export
setMethod("getPiRNAID", "piRNAGene", function(object) object@piRNAID)
