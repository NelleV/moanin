# Normalization utility functions
#

#' Utility function to filter out low-expressed genes
#'
#' @param counts n by p matrix containing the count data.
#' @param min_counts integer, corresponding to the minimum number of counts for
#'  the gene to be considered expressed in a sample. Default is 20.
#' @param min_samples integer, corresponding to the minimum of samples for a
#'  gene to be expressed to be included in downstream analyses.
#' @return The filtered counts matrix
expression_filtering = function(counts, min_counts=20, min_samples=3){
    rows_to_keep = apply(counts, 1, function(x){
            sum(x > min_counts) > min_samples
    })
    counts = counts[rows_to_keep, ]
    return(counts)
}

