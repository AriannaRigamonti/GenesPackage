% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_getDescription.R
\name{getDescription}
\alias{getDescription}
\alias{getDescription,Gene-method}
\title{Get the gene description}
\usage{
getDescription(object)

\S4method{getDescription}{Gene}(object)
}
\arguments{
\item{object}{A Gene object.}
}
\value{
The gene description.
}
\description{
This function retrieves the description of a gene.
}
\examples{
gene <- createGene(1L, "SYMBOL", "Gene Name", "Description", "chr1", 1, 1000, "+", list())
getDescription(gene)
}
