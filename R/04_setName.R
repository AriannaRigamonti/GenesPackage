#' Set the gene name
#'
#' This method sets the name of a Gene object.
#'
#' @param object A Gene object.
#' @param value The new gene name.
#' @return The updated Gene object.
#' @title Set the gene name
#' @name setName
#' @aliases setName<- setName<-,Gene-method
#' @rdname setName
#' @export
#' @examples
#' gene <- createGene(
#'   1L, "SYMBOL", "Gene Name", "Description",
#'   "chr1", 1, 1000, "+", list()
#' )
#' setName(gene) <- "New Gene Name"
#' getName(gene)
setGeneric("setName<-", function(object, value) standardGeneric("setName<-"))

#' @rdname setName
#' @export
setMethod("setName<-", "Gene", function(object, value) {
  object@name <- value
  object
})
