#' Set the gene description
#'
#' This function sets the description of a gene.
#'
#' @param object A Gene object.
#' @param value The new gene description.
#' @return The updated Gene object.
#' @title Set the gene description
#' @name setDescription
#' @aliases setDescription<- setDescription<-,Gene-method
#' @rdname setDescription
#' @export
#' @examples
#' gene <- createGene(1L, "SYMBOL", "Gene Name", "Description", "chr1", 1, 1000, "+", list())
#' setDescription(gene) <- "New Description"
#' getDescription(gene)
setGeneric("setDescription<-", function(object, value) standardGeneric("setDescription<-"))

#' @rdname setDescription
#' @export
setMethod("setDescription<-", "Gene", function(object, value) {
  object@description <- value
  object
})
