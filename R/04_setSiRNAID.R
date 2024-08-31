#' Set the siRNA ID
#' @name setSiRNAID<-
#' @title Set the siRNA ID
#' @param object A siRNAGene object.
#' @param value The new siRNA ID.
#' @return The updated siRNAGene object.
#' @export
#' @aliases setSiRNAID setSiRNAID<-,siRNAGene-method
#' @examples
#' sirna_gene <- createSiRNAGene(
#'   4L, "SYMBOL_SI", "siRNA Name",
#'   "siRNA Description", "chr1", 1, 1000, "+",
#'   list(), "sirna1", "SIRNA_SEQ"
#' )
#' setSiRNAID(sirna_gene) <- "new_sirna1"
#' getSiRNAID(sirna_gene)
setGeneric("setSiRNAID<-", function(object, value) standardGeneric("setSiRNAID<-"))
setMethod("setSiRNAID<-", "siRNAGene", function(object, value) {
  object@siRNAID <- value
  object
})
