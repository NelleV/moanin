#' Small data set for running examples
#'
#' @name exampleData
#' @docType data
#' @format Two objects are loaded, a data frame of expression of 500 genes by 84
#'   samples (\code{testData}), and a data frame with meta information on those
#'   84 samples (\code{testMeta}).
#' @details This data is a subset of the full time course data available as
#'   \code{\link[timecoursedata]{shoemaker2015}} and is only provided for the
#'   purpose of running examples, and not for biological meaning. Users should
#'   refer to the full data set.
#' @keywords data
#' @examples
#' #code used to create data:
#' \dontrun{
#' library(timecoursedata)
#' data(shoemaker2015)
#' testData<-shoemaker2015$data[1:500,]
#' whSamples<-which(shoemaker2015$meta$Group %in% c("C","K"))
#' testData<-testData[,whSamples]
#' testMeta<-shoemaker2015$meta[whSamples,]
#' save(list=c("testData","testMeta"),file="data/exampleData.rda")
#' }
NULL

