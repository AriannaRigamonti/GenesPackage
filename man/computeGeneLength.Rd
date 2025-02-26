% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_computeGeneLength.R
\name{computeGeneLength}
\alias{computeGeneLength}
\alias{computeGeneLength,Gene-method}
\title{Compute the length of the gene}
\usage{
computeGeneLength(object)

\S4method{computeGeneLength}{Gene}(object)
}
\arguments{
\item{object}{A Gene object.}
}
\value{
The length of the gene as an integer.
}
\description{
This function calculates the length of a gene based on its genomic structure.
The length is computed as the sum of the widths of the ranges in the gene's structure.
}
\examples{
gene <- createGene(
  ID = 1L, symbol = "SYMBOL", name = "Gene Name",
  description = "Description", chr = "chr1", start = 1,
  end = 1000, strand = "+", product = list()
)
computeGeneLength(gene)
}
