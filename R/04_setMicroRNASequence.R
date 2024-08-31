#' Set the microRNA sequence
#'
#' This method sets the microRNA sequence of a microRNAGene object.
#'
#' @param object A microRNAGene object.
#' @param value The new microRNA sequence.
#' @return The updated microRNAGene object.
#' @title Set the microRNA sequence
#' @name setMicroRNASequence
#' @aliases setMicroRNASequence<- setMicroRNASequence<-,microRNAGene-method
#' @rdname setMicroRNASequence
#' @export
#' @examples
#' mirna_gene <- createMicroRNAGene(
#'   3L, "SYMBOL_MIR", "MicroRNA Name",
#'   "MicroRNA Description", "chr1", 1, 1000, "+",
#'   list(), "mirna1", "SEED_SEQ"
#' )
#' setMicroRNASequence(mirna_gene) <- "NEW_SEED_SEQ"
#' getMicroRNASequence(mirna_gene)
setGeneric("setMicroRNASequence<-", function(object, value) standardGeneric("setMicroRNASequence<-"))

#' @rdname setMicroRNASequence
#' @export
setMethod("setMicroRNASequence<-", "microRNAGene", function(object, value) {
  object@microRNASequence <- value
  object
})
