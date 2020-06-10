% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splines.R
\name{pvalues_fisher_method}
\alias{pvalues_fisher_method}
\title{Fisher's method to combine pvalues}
\usage{
pvalues_fisher_method(pvalues)
}
\arguments{
\item{pvalues}{a matrix of pvalues, with columns corresponding to different
tests or sources of p-values, and rows corresponding to the genes from
which the p-values come.}
}
\value{
a vector of p-values, one for each row of \code{pvalues}, that is the
  result of Fisher's combined probability test applied to the p-values in
  that row.
}
\description{
Combines all p-value per rows.
}