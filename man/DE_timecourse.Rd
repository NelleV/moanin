% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DE_timecourse.R
\name{DE_timecourse}
\alias{DE_timecourse}
\alias{DE_timecourse,Moanin-method}
\title{Run spline models and test for DE of contrasts.}
\usage{
\S4method{DE_timecourse}{Moanin}(
  object,
  contrasts,
  center = FALSE,
  statistic = c("ftest", "lrt"),
  use_voom_weights = TRUE
)
}
\arguments{
\item{object}{An object of class \code{\link{Moanin}}, an object containing
all related information for time course data and the splines model that
will be used (if applicable). See \code{\link{create_moanin_model}} for
more details.}

\item{contrasts}{Contrasts, either provided as a vector of strings, or a
matrix of contrasts coefficients obtained using
\code{\link[limma]{makeContrasts}} from the package \code{limma}. If given
as a character string, will be passed to \code{\link[limma]{makeContrasts}}
to be converted into such a matrix.}

\item{center}{boolean, whether to center the data matrix}

\item{statistic}{Which test statistic to use, a likelihood ratio statistic or
a F-test.}

\item{use_voom_weights}{boolean, optional, default: TRUE. 
Whether to use voom weights. See details.}
}
\value{
A \code{data.frame} with two columns for each of the contrasts given
  in \code{contrasts}, corresponding to the raw p-value of the contrast for
  that gene (\code{_pval}) and the adjusted p-value (\code{_qval}). The
  adjusted p-values are FDR-adjusted based on the Benjamini-Hochberg method,
  as implemented in \code{\link[stats]{p.adjust}}. The adjustment is done
  across all p-values for all contrasts calculated.
}
\description{
Run spline models and test for DE of contrasts.
}
\details{
The implementation of the spline fit and the calculation of p-values
  was based on code from \code{\link[edge]{edge}}, and expanded to enable
  handling of comparisons of groups via contrasts. The code assumes that the \code{Moanin} object was created via either a formula or a basis where a different spline was fit for each \code{group_variable} and thus the contrasts are comparisons of those spline fits. If the \code{Moanin} object was created via user-provided basis matrix or formula, then the user should take a great deal of caution in using this code, as the degrees of freedom for the tests of significance cannot be verified to be correct.

If \code{use_voom_weights=TRUE}, then before fitting splines to each gene,
voom weights are calculated from \code{assay(object)}:
\preformatted{
   y <- edgeR::DGEList(counts=assay(object))
   y <- edgeR::calcNormFactors(y, method="upperquartile")
   v <- limma::voom(y, design, plot=FALSE)
   weights <- v$weights
}
The design matrix for the voom weights is based on the formula
\code{~Group + Timepoint +0}
where Group and Timepoint are replaced with the user-defined values where appropriate. 
These weights are given to the \code{lm.fit} which fits the spline coefficients.
This workflow assumes that the input to the \code{Moanin} object were counts.

If the user set \code{log_transform=TRUE} in the creation of the
  \code{Moanin} object, the splines will be fit to the log of the input data,
  and not directly to the input data. This is independent of whether the user
  chooses \code{use_voom_weights}.
}
\examples{
data(exampleData)
moanin <- create_moanin_model(data=testData, meta=testMeta)
deTimecourse=DE_timecourse(moanin, 
   contrasts="K-C", use_voom_weights=FALSE)
head(deTimecourse)
}
\seealso{
\code{\link[limma]{makeContrasts}}, \code{\link{create_moanin_model}},
  \code{\link{DE_timepoints}}, \code{\link[edge]{edge}}
}
