#' Get the gene ID
#'
#' This function retrieves the ID of a gene.
#'
#' @param object A Gene object.
#' @return The gene ID.
#' @export
#' @aliases getID getID,Gene-method
#' @examples
#' gene <- createGene(1L, "SYMBOL", "Gene Name", "Description", "chr1", 1, 1000, "+", list())
#' getID(gene)
setGeneric("getID", function(object) standardGeneric("getID"))

#' @rdname getID
#' @export
setMethod("getID", "Gene", function(object) object@ID)
