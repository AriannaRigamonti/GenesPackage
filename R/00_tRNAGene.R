#' tRNAGene class
#'
#' A class representing a tRNA gene.
#'
#' @slot tRNAID character. The tRNA ID.
#' @slot tRNASequence character. The tRNA sequence.
#' @exportClass tRNAGene
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom IRanges IRanges
#' @importFrom methods new
#' @examples
#' gr <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges::IRanges(start = 1, end = 1000)
#' )
#' trna_gene <- new("tRNAGene",
#'   ID = 6L, symbol = "SYMBOL_T",
#'   name = "tRNA Name", description = "tRNA Description",
#'   structure = gr, product = list(), tRNAID = "trna1",
#'   tRNASequence = "TRNA_SEQ"
#' )
setClass("tRNAGene",
  contains = "Gene",
  slots = list(
    tRNAID = "character",
    tRNASequence = "character"
  )
)
