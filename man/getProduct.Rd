% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_getProduct.R
\name{getProduct}
\alias{getProduct}
\alias{getProduct,Gene-method}
\title{Get the gene product}
\usage{
getProduct(object)

\S4method{getProduct}{Gene}(object)
}
\arguments{
\item{object}{A Gene object.}
}
\value{
The gene product.
}
\description{
This function retrieves the product of a gene.
}
\examples{
gene <- createGene(1L, "SYMBOL", "Gene Name", "Description", "chr1", 1, 1000, "+", list())
getProduct(gene)
}
