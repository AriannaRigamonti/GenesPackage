% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/05_summary-snoRNAGene.R
\name{summary-snoRNAGene}
\alias{summary-snoRNAGene}
\alias{summary,snoRNAGene-method}
\title{Summary Method for snoRNAGene Class}
\usage{
\S4method{summary}{snoRNAGene}(object)
}
\arguments{
\item{object}{A snoRNAGene object.}
}
\value{
A summary of the snoRNA gene information.
}
\description{
Summary method for snoRNAGene objects
}
\details{
Provides a summary of the snoRNA gene information.
}
\examples{
snorna_gene <- createSnoRNAGene(
  5L, "SYMBOL_SNO", "snoRNA Name",
  "snoRNA Description", "chr1", 1, 1000,
  "+", list(), "snorna1", "SNORNA_SEQ"
)
summary(snorna_gene)
}
