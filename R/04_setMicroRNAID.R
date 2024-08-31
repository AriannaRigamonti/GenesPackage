#' Set the microRNA ID
#'
#' This method sets the microRNA ID of a microRNAGene object.
#'
#' @param object A microRNAGene object.
#' @param value The new microRNA ID.
#' @return The updated microRNAGene object.
#' @title Set the microRNA ID
#' @name setMicroRNAID
#' @aliases setMicroRNAID<- setMicroRNAID<-,microRNAGene-method
#' @rdname setMicroRNAID
#' @export
#' @examples
#' mirna_gene <- createMicroRNAGene(
#'   3L, "SYMBOL_MIR", "MicroRNA Name",
#'   "MicroRNA Description", "chr1", 1, 1000,
#'   "+", list(), "mirna1", "SEED_SEQ"
#' )
#' setMicroRNAID(mirna_gene) <- "new_mirna1"
#' getMicroRNAID(mirna_gene)
setGeneric("setMicroRNAID<-", function(object, value) standardGeneric("setMicroRNAID<-"))

#' @rdname setMicroRNAID
#' @export
setMethod("setMicroRNAID<-", "microRNAGene", function(object, value) {
  object@microRNAID <- value
  object
})
