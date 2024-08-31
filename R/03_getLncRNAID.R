#' Get the lncRNA ID
#'
#' This function retrieves the ID of a lncRNA gene.
#'
#' @param object A lncRNAGene object.
#' @return The lncRNA ID.
#' @export
#' @aliases getLncRNAID getLncRNAID,lncRNAGene-method
#' @examples
#' lncrna_gene <- createLncRNAGene(
#'   2L, "SYMBOL_LNC", "LncRNA Name",
#'   "LncRNA Description", "chr1", 1, 1000, "+",
#'   list(), "lncrna1", "RNA_SEQUENCE"
#' )
#' getLncRNAID(lncrna_gene)
setGeneric("getLncRNAID", function(object) standardGeneric("getLncRNAID"))

#' @rdname getLncRNAID
#' @export
setMethod("getLncRNAID", "lncRNAGene", function(object) object@lncRNAID)
