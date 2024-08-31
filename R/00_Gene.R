#' Define the Gene class
#' Gene class
#'
#' A class representing a generic gene.
#'
#' @slot ID integer The gene ID (e.g., Ensembl ID or NCBI gene ID).
#' @slot symbol character. The gene symbol.
#' @slot name character. The gene name.
#' @slot description character. The gene description.
#' @slot structure GRanges. The gene structure.
#' @slot product list. The gene product.
#' @exportClass Gene
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom IRanges IRanges
#' @import methods
#' @examples
#' gr <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges::IRanges(start = 1, end = 1000)
#' )
#' gene <- new("Gene",
#'   ID = 1L,
#'   symbol = "SYMBOL",
#'   name = "Gene Name",
#'   description = "Description",
#'   structure = gr,
#'   product = list()
#' )
setClass("Gene",
  slots = list(
    ID = "integer",
    symbol = "character",
    name = "character",
    description = "character",
    structure = "GRanges",
    product = "list"
  )
)
