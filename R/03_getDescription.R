#' Get the gene description
#'
#' This function retrieves the description of a gene.
#'
#' @param object A Gene object.
#' @return The gene description.
#' @export
#' @aliases getDescription getDescription,Gene-method
#' @examples
#' gene <- createGene(1L, "SYMBOL", "Gene Name", "Description", "chr1", 1, 1000, "+", list())
#' getDescription(gene)
setGeneric("getDescription", function(object) standardGeneric("getDescription"))

#' @rdname getDescription
#' @export
setMethod("getDescription", "Gene", function(object) object@description)
