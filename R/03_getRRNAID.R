#' Get the rRNA ID
#'
#' This function retrieves the ID of a rRNA gene.
#'
#' @param object A rRNAGene object.
#' @return The rRNA ID.
#' @export
#' @aliases getRRNAID getRRNAID,rRNAGene-method
#' @examples
#' rrna_gene <- createRRNAGene(
#'   7L, "SYMBOL_R", "rRNA Name",
#'   "rRNA Description", "chr1", 1, 1000, "+", list(),
#'   "rrna1", "RRNA_SEQ"
#' )
#' getRRNAID(rrna_gene)
setGeneric("getRRNAID", function(object) standardGeneric("getRRNAID"))

#' @rdname getRRNAID
#' @export
setMethod("getRRNAID", "rRNAGene", function(object) object@rRNAID)
