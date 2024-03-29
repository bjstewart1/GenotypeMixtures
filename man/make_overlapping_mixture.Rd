% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gt_mix_functions.R
\name{make_overlapping_mixture}
\alias{make_overlapping_mixture}
\title{make_overlapping_mixture}
\usage{
make_overlapping_mixture(n_mixtures, n_genotypes, density = 1)
}
\arguments{
\item{n_mixtures}{the number of mixtures to include}

\item{n_genotypes}{the number of genotypes to include}

\item{density}{the density of the design - numeric value from 0->1}
}
\value{
matrix of mixtures by genotypes
}
\description{
Function produces an experimental design with varying degrees of sparsity.
}
\examples{
\donttest{
make_overlapping_mixture()
}
}
