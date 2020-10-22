setGeneric("estimate_log_fold_change", 
           function(object,...) { standardGeneric("estimate_log_fold_change")})

#' Estimates log fold change
#'
#' @inheritParams DE_timecourse
#' @param contrasts The contrasts to consider
#' @param method method for calculating the log-fold change. See details. 
#' @details The following methods exist for calculating the log-fold change
#'   between conditions over time (default is "timecourse"):
#' \itemize{
#' \item{\code{timely}}{The log-fold change for each individual timepoint
#' (\eqn{lfc(t)})}
#' \item{\code{timecourse}}{The average absolute per-week fold-change,
#' multiplied by the sign of the average per-week fold-change.}
#' \item{\code{sum}}{Sum of per-week log fold change, over all timepoints}
#' \item{\code{max}}{Max of per-week log fold change, over all timepoints}
#' \item{\code{abs_sum}}{Sum of the absolute value of the per-week log fold
#' change, over all timepoints}
#' \item{\code{abs_squared_sum}}{Sum of the square value of the per-week log
#' fold change, over all timepoint}
#' \item{\code{min}}{Min of per-week log fold change, over all timepoints}
#' }
#' @details If the user set \code{log_transform=TRUE} in the creation of the
#'   \code{Moanin} object, the data will be log transformed before calculating
#'   the fold-change.
#' @return A data.frame giving the estimated log-fold change for each gene
#'   (row). For all methods except for "timely", the data frame will consist of
#'   one column for each value of the argument \code{contrasts}. For "timely"
#'   there will be one column for each timepoint and contrast combination.
#' @name estimate_log_fold_change
#' @aliases estimate_log_fold_change,Moanin-method
#' @examples 
#' data(exampleData)
#' moanin <- create_moanin_model(data=testData,meta=testMeta)
#' estsTimely <- estimate_log_fold_change(moanin, 
#'    contrasts=c("K-C"), method="timely")
#' head(estsTimely)
#' estsTimecourse <- estimate_log_fold_change(moanin, 
#'    contrasts=c("K-C"),method="timecourse")
#' head(estsTimecourse)
#' @export
setMethod("estimate_log_fold_change", "Moanin",
          function(object, 
        contrasts, 
        method=c("timecourse", "sum", "max", "timely", 
                 "abs_sum", "abs_squared_sum", "min")){

    tpVar <- time_variable_name(object)
    gpVar <- group_variable_name(object)
    method <- match.arg(method)
    contrasts <- is_contrasts(contrasts, object)
    
    if(method == "timely"){
        log_fold_changes <- lfc_per_time(object, contrasts)
    }else if(method %in% c("sum", "max", "min", "abs_sum", "abs_squared_sum",
                           "timecourse", "epicon")){
        timely_lfc <- lfc_per_time(object, contrasts)
        timely_lfc_meta <- reconstruct_meta_from_lfc(timely_lfc, split_char=":",
                    group_variable_name=gpVar,time_variable_name=tpVar)
        log_fold_changes <- data.frame(row.names=row.names(object))
        for(contrast in colnames(contrasts)){
            mask <- (timely_lfc_meta[,gpVar] == contrast) & 
                !is.na(matrixStats::colSums2(timely_lfc))
            if(method == "max"){
                log_fold_changes[, contrast] <- matrixStats::rowMaxs(abs(timely_lfc[, mask]))
            } else if(method == "min") {
                log_fold_changes[, contrast] <- matrixStats::rowMins(abs(timely_lfc[, mask])) 
            }else if(method == "abs_sum"){
                log_fold_changes[, contrast] <- matrixStats::rowSums2(abs(timely_lfc[, mask]))
            }else if(method == "abs_squared_sum"){
                log_fold_changes[, contrast] <- matrixStats::rowSums2(timely_lfc[, mask]**2)
            }else if(method == "epicon" || method == "timecourse"){
                log_fold_changes[, contrast] <- (
                    matrixStats::rowMeans2(abs(timely_lfc[, mask])) * 
                        sign(matrixStats::rowSums2(timely_lfc[, mask])))
            }else if(method == "sum"){
                log_fold_changes[, contrast] <- matrixStats::rowSums2(timely_lfc[, mask])
            }
        }
        
    }
    
    return(log_fold_changes)
    
}
)

# XXX helper function to reconstruct metadat from adat
reconstruct_meta_from_lfc <- function(data_per_time, group_variable_name,
                                     time_variable_name,split_char="."){
    meta_per_time <- t(
        as.data.frame(strsplit(colnames(data_per_time), split_char, 
                               fixed=TRUE)))
    row.names(meta_per_time) <- colnames(data_per_time)
    colnames(meta_per_time) <- c(group_variable_name, time_variable_name)
    meta_per_time <- as.data.frame(meta_per_time)
    meta_per_time[, time_variable_name] <- 
        as.numeric(meta_per_time[, time_variable_name])
    return(meta_per_time)
}


lfc_per_time <- function(object, contrasts){
    tpVar <- time_variable_name(object)
    gpVar <- group_variable_name(object)
    group_variable(object) <- as.factor(group_variable(object))
    averaged_data <- average_replicates(object)
    averaged_meta <- reconstruct_meta_from_lfc(averaged_data, split_char=":",
                                group_variable_name=gpVar,
                                time_variable_name=tpVar)
    
    averaged_meta[,tpVar] <- as.factor(averaged_meta[,tpVar])
    sample_coefficients <- vapply(
        averaged_meta[,gpVar], 
        function(x) return(contrasts[x, ]),
        numeric(ncol(contrasts)))

    if(is.null(dim(sample_coefficients))){
        sample_coefficients <- as.matrix(sample_coefficients)
    }else{
        sample_coefficients <- t(sample_coefficients)
    }

    colnames(sample_coefficients) <- colnames(contrasts)
    
    log_fold_changes <- matrix(NA, dim(object)[1],
                     dim(contrasts)[2]*length(unique(time_variable(object))))
    row.names(log_fold_changes) <- row.names(object)
    colnames(log_fold_changes) <- vapply(
        colnames(sample_coefficients),
        FUN=function(x){vapply(unique(averaged_meta[,tpVar]), 
                               FUN=function(t){paste0(x, ":", t)},
                               character(1))},
        character(length(unique(averaged_meta[, tpVar]))))
    
    for(column in colnames(sample_coefficients)){
        
        sample_coefficient <- as.vector(unlist(sample_coefficients[, column]))
        coef_data <- t(t(averaged_data) * sample_coefficient)
        for(timepoint in unique(averaged_meta[,tpVar])){
            mask <- averaged_meta[,tpVar] == timepoint
            if(sum(mask) != dim(contrasts)[1]){
                next
            }
            colname <- paste0(column, ":", as.character(timepoint))
            log_fold_changes[, colname] <- matrixStats::rowSums2(coef_data[, mask])
        }
    }
    return(log_fold_changes)
}


average_replicates <- function(object){
    timepoint_group <- 
        droplevels(group_variable(object):as.factor(time_variable(object)))
    all_levels <- levels(timepoint_group)
    
    # selecting certain columns sometimes returns vectors and sometimes 
    # matrices depending on whether the user wants a single column or several
    # columns. rowMeans requires a matrix. Thus, here's a small implementation
    # that will work in all cases.
    .row_means <- function(data){
        if(is.vector(data)){
            return(data)
        }else{
            return(matrixStats::rowMeans2(data))
        }
    }
    y <- get_log_data(object)
    if(inherits(y, "DataFrame")){
        y <- data.matrix(y)
    }
    replicate_averaged <- vapply(
        all_levels,
        function(m){.row_means(y[, timepoint_group==m])},
        numeric(nrow(y)))
    return(replicate_averaged)
}

    
