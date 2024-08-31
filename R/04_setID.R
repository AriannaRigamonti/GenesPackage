#' Set the gene ID
#'
#' This method sets the ID of a Gene object.
#'
#' @param object A Gene object.
#' @param value The new gene ID.
#' @return The updated Gene object.
#' @title Set the gene ID
#' @name setID
#' @aliases setID<- setID<-,Gene-method
#' @rdname setID
#' @export
#' @examples
#' gene <- createGene(1L, "SYMBOL", "Gene Name", "Description", "chr1", 1, 1000, "+", list())
#' setID(gene) <- 2L
#' getID(gene)
setGeneric("setID<-", function(object, value) standardGeneric("setID<-"))

#' @rdname setID
#' @export
setMethod("setID<-", "Gene", function(object, value) {
  object@ID <- value
  object
})
