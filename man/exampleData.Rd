% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exampleData.R
\docType{data}
\name{exampleData}
\alias{exampleData}
\alias{testData}
\alias{testMeta}
\alias{testGenesGO}
\title{Small data set for running examples}
\format{
Three objects are loaded, a data frame of expression of 500 genes by
  84 samples (\code{testData}), a data frame with meta information on those
  84 samples (\code{testMeta}), and a data frame giving the GOID of the genes
  in \code{testData}.
}
\description{
Small data set for running examples
}
\details{
This data is a subset of the full time course data available as
  \code{\link[timecoursedata]{shoemaker2015}} and is only provided for the
  purpose of running examples, and not for biological meaning. Users should
  refer to the full data set.

The rownames of \code{testData} are RefSeq.
}
\examples{
#code used to create data:
\dontrun{
library(timecoursedata)
data(shoemaker2015)
testData<-shoemaker2015$data[1:500,]
whSamples<-which(shoemaker2015$meta$Group \%in\% c("C","K","M"))
testData<-testData[,whSamples]
testMeta<-droplevels(shoemaker2015$meta[whSamples,])
library(biomaRt)
ensembl = useMart("ensembl")
ensembl = useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl", mart=ensembl)
testGenesGO = getBM(attributes=c("go_id", "refseq_mrna"),
              values=rownames(testData),
              filters="refseq_mrna",
              mart=ensembl)
save(list=c("testData","testMeta","testGenesGO"),file="data/exampleData.rda")
}
}
\keyword{data}
