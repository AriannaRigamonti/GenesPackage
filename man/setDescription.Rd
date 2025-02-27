% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_setDescription.R
\name{setDescription}
\alias{setDescription}
\alias{setDescription<-}
\alias{setDescription<-,Gene-method}
\title{Set the gene description}
\usage{
setDescription(object) <- value

\S4method{setDescription}{Gene}(object) <- value
}
\arguments{
\item{object}{A Gene object.}

\item{value}{The new gene description.}
}
\value{
The updated Gene object.
}
\description{
Set the gene description
}
\details{
This function sets the description of a gene.
}
\examples{
gene <- createGene(1L, "SYMBOL", "Gene Name", "Description", "chr1", 1, 1000, "+", list())
setDescription(gene) <- "New Description"
getDescription(gene)
}
