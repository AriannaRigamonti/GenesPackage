% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_getPiRNAID.R
\name{getPiRNAID}
\alias{getPiRNAID}
\alias{getPiRNAID,piRNAGene-method}
\title{Get the piRNA ID}
\usage{
getPiRNAID(object)

\S4method{getPiRNAID}{piRNAGene}(object)
}
\arguments{
\item{object}{A piRNAGene object.}
}
\value{
The piRNA ID.
}
\description{
This function retrieves the ID of a piRNA gene.
}
\examples{
pirna_gene <- createPiRNAGene(
  8L, "SYMBOL_PI", "piRNA Name",
  "piRNA Description", "chr1", 1, 1000, "+",
  list(), "pirna1", "PIRNA_SEQ"
)
getPiRNAID(pirna_gene)
}
