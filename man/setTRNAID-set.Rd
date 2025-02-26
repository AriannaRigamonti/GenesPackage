% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_setTRNAID.R
\name{setTRNAID<-}
\alias{setTRNAID<-}
\alias{setTRNAID}
\alias{setTRNAID<-,tRNAGene-method}
\title{Set the tRNA ID}
\usage{
setTRNAID(object) <- value
}
\arguments{
\item{object}{A tRNAGene object.}

\item{value}{The new tRNA ID.}
}
\value{
The updated tRNAGene object.
}
\description{
Set the tRNA ID
}
\examples{
trna_gene <- createTRNAGene(
  6L, "SYMBOL_T", "tRNA Name", "tRNA
                            Description", "chr1", 1, 1000, "+", list(),
  "trna1", "TRNA_SEQ"
)
setTRNAID(trna_gene) <- "new_trna1"
getTRNAID(trna_gene)
}
