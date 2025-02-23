% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_setSnoRNASequence.R
\name{setSnoRNASequence<-}
\alias{setSnoRNASequence<-}
\alias{setSnoRNASequence}
\alias{setSnoRNASequence<-,snoRNAGene-method}
\title{Set the snoRNA sequence}
\usage{
setSnoRNASequence(object) <- value
}
\arguments{
\item{object}{A snoRNAGene object.}

\item{value}{The new snoRNA sequence.}
}
\value{
The updated snoRNAGene object.
}
\description{
Set the snoRNA sequence
}
\examples{
snorna_gene <- createSnoRNAGene(
  5L, "SYMBOL_SNO", "snoRNA Name",
  "snoRNA Description", "chr1", 1, 1000, "+",
  list(), "snorna1", "SNORNA_SEQ"
)
setSnoRNASequence(snorna_gene) <- "NEW_SNORNA_SEQ"
getSnoRNASequence(snorna_gene)
}
