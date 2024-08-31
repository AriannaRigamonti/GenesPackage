#' Set the piRNA ID
#'
#' This method sets the piRNA ID of a piRNAGene object.
#'
#' @param object A piRNAGene object.
#' @param value The new piRNA ID.
#' @return The updated piRNAGene object.
#' @title Set the piRNA ID
#' @name setPiRNAID
#' @aliases setPiRNAID<- setPiRNAID<-,piRNAGene-method
#' @rdname setPiRNAID
#' @export
#' @examples
#' pirna_gene <- createPiRNAGene(
#'   8L, "SYMBOL_PI", "piRNA Name",
#'   "piRNA Description", "chr1", 1, 1000,
#'   "+", list(), "pirna1", "PIRNA_SEQ"
#' )
#' setPiRNAID(pirna_gene) <- "new_pirna1"
#' getPiRNAID(pirna_gene)
setGeneric("setPiRNAID<-", function(object, value) {
  standardGeneric("setPiRNAID<-")
})

#' @rdname setPiRNAID
#' @export
setMethod("setPiRNAID<-", "piRNAGene", function(object, value) {
  object@piRNAID <- value
  object
})
