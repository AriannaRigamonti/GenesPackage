#' Set the gene product
#'
#' This method sets the product of a Gene object.
#'
#' @param object A Gene object.
#' @param value The new gene product.
#' @return The updated Gene object.
#' @title Set the gene product
#' @name setProduct
#' @aliases setProduct<- setProduct<-,Gene-method
#' @rdname setProduct
#' @export
#' @examples
#' gene <- createGene(
#'   1L, "SYMBOL", "Gene Name", "Description",
#'   "chr1", 1, 1000, "+", list()
#' )
#' setProduct(gene) <- list("New Product")
#' getProduct(gene)
setGeneric("setProduct<-", function(object, value) standardGeneric("setProduct<-"))

#' @rdname setProduct
#' @export
setMethod("setProduct<-", "Gene", function(object, value) {
  object@product <- value
  object
})
