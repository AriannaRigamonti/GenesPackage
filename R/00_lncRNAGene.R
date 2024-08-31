#' lncRNAGene class
#'
#' A class representing a long non-coding RNA gene.
#'
#' @slot lncRNAID character The lncRNA ID.
#' @slot RNASequence character. The RNA sequence.
#' @exportClass lncRNAGene
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom IRanges IRanges
#' @importFrom methods new
#' @examples
#' gr <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges::IRanges(start = 1, end = 1000)
#' )
#' lncrna_gene <- new("lncRNAGene",
#'   ID = 2L, symbol = "SYMBOL_LNC",
#'   name = "LncRNA Name", description = "LncRNA Description",
#'   structure = gr, product = list(), lncRNAID = "lncrna1",
#'   RNASequence = "RNA_SEQUENCE"
#' )
setClass("lncRNAGene",
  contains = "Gene",
  slots = list(
    lncRNAID = "character",
    RNASequence = "character"
  )
)
