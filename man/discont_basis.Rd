% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/discont_basis.R
\name{discont_basis}
\alias{discont_basis}
\title{Provides set of basis functions on either side of a time point, allowing for a discontinuity in the fitted functions}
\usage{
discont_basis(
  timepoints,
  discont_point,
  knots = NULL,
  dfPre = NULL,
  dfPost = dfPre,
  degree = 3,
  intercept = TRUE,
  type = c("ns", "bs")
)
}
\arguments{
\item{timepoints}{vector of numeric timepoints for which the splines basis will be evaluated}

\item{discont_point}{a single numeric value that represents where the discontinuity should be}

\item{knots}{passed to \code{ns} or \code{bs}. If not NULL, should give knots on either side of \code{discon_point} as single vector -- they will be separated in the call to \code{discon_point}}

\item{dfPre}{the df for the basis functions defined before the discontinuity point}

\item{dfPost}{the df for the basis functions defined after the discontinuity point}

\item{degree}{passed to \code{bs} (if applicable)}

\item{intercept}{Whether to include an intercept (vector of all 1s) for each side of the discontinuity. Note this is different than the argument \code{intercept} of either \code{bs} or \code{ns}, which is set to \code{FALSE}.}

\item{type}{either "ns" or "bs" indicating which splines basis function should be used.}
}
\description{
Provides set of basis functions on either side of a time point, allowing for a discontinuity in the fitted functions
}
\examples{
x<-seq(0,10,length=100)
basis<-discont_basis(x,discont_point=3, dfPre=3, dfPost=4, intercept=TRUE)
# Plot of the basis functions
par(mfrow=c(3,3))
for(i in 1:ncol(basis)){
   plot(x,basis[,i],type="l")
   abline(v=3,lty=2)
}
# Use it in a moanin_model object instead of ns/bs:
data(exampleData)
moanin <- create_moanin_model(data=testData, meta=testMeta,
    spline_formula=~Group:discont_basis(Timepoint,dfPre=3,
        dfPost=3,discont=20,intercept=TRUE)+0,
    degrees_of_freedom=6)
}
