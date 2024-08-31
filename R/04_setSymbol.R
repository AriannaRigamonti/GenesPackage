#' Set the gene symbol
#'
#' This method sets the symbol of a Gene object.
#'
#' @param object A Gene object.
#' @param value The new gene symbol.
#' @return The updated Gene object.
#' @title Set the gene symbol
#' @name setSymbol
#' @aliases setSymbol<- setSymbol<-,Gene-method
#' @rdname setSymbol
#' @export
#' @examples
#' gene <- createGene(
#'   1L, "SYMBOL", "Gene Name", "Description",
#'   "chr1", 1, 1000, "+", list()
#' )
#' setSymbol(gene) <- "NewSymbol"
#' getSymbol(gene)
setGeneric("setSymbol<-", function(object, value) standardGeneric("setSymbol<-"))

#' @rdname setSymbol
#' @export
setMethod("setSymbol<-", "Gene", function(object, value) {
  object@symbol <- value
  object
})
