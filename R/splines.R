setGeneric("rescale_values", function(object,...) { 
    standardGeneric("rescale_values")})

#' Fit splines to each gene of data matrix, used by DE_timecourse
#' @param design_matrix design matrix (containing evaluated basis matrix, but potentially other variables) to fit
#' @param data a matrix of data to fix splines to
#' @param weights A matrix of weights, of the same dimension as \code{data}.
#' @details Needed to allow for fitting weights, in which case for loop over every row/gene. Otherwise, just a call to lm.fit.
#' @return matrix of the coefficients for each basis function, each row of the
#'   matrix containing the coefficients for the corresponding gene in
#'   \code{data}.
#' @keywords internal
fit_splines <- function(design_matrix, data, weights=NULL){
    n <- ncol(design_matrix)
    nr <- nrow(data)
    if(inherits(data,"DataFrame")) data<-data.matrix(data)
    
    if(!is.null(weights)){
        beta <- matrix(nrow=nr, ncol=n)
        ##FIXME: can't we multiply data by weights and then make call to lm.fit -- would be faster than looping over all genes, no?
        for(i in seq_len(nr)){
            beta[i,] <- stats::lm.wfit(design_matrix, data[i,], weights[i,])$coefficients
        }
        row.names(beta) <- row.names(data)
    }else{
        beta <- t(stats::lm.fit(design_matrix, t(data))$coefficients)
    }
    return(beta)
}

#' Get fitted values for splines for each gene
#'
#' @inheritParams DE_timecourse
#' @param meta_prediction 
#' @param data a matrix with data (required doesn't pull from moanin_model)
#' @return a matrix of the fitted y values, with dimensions the same as
#'   \code{data}
#'
#' @keywords internal
fit_predict_splines <- function(data, moanin_model, 
                                meta_prediction=NULL){
    basis <- basis_matrix(moanin_model)
    gpVar <- group_variable_name(moanin_model)
    tpVar <- time_variable_name(moanin_model)

    if(inherits(data,"DataFrame")){
        data <- data.matrix(data)
    }

    if(is.null(meta_prediction)){
        y_fitted <- t(stats::lm.fit(basis, t(data))$fitted.values)
    }else{
        degrees_of_freedom <- degrees_of_freedom(moanin_model)
        fitting_data <- t(as.matrix(data))
        formula_data <- list(
            "Group"=group_variable(moanin_model),
            "Timepoint"=time_variable(moanin_model),
            "fitting_data"=fitting_data,
            "degrees_of_freedom"=degrees_of_freedom(moanin_model))

        names(formula_data)[c(1, 2)] <- c(gpVar, tpVar)
        updated_formula <- stats::update(spline_formula(moanin_model), 
                                        fitting_data ~ .)
        model <- stats::lm(updated_formula, formula_data)
        y_fitted <- stats::predict(model, meta_prediction)
    }
    return(y_fitted)
}


#' Create prediction meta data from splines model
#'
#' @inheritParams DE_timecourse
#' @param num_timepoints integer, optional, default: 100.
#'  Number of timepoints to use for the prediction metadata
#'
#' @return moanin_model
#' @keywords internal
create_meta_prediction <- function(moanin_model, num_timepoints=100){
    # Create moanin_model for prediction
    timepoints_pred <- NULL
    groups_pred <- NULL
    gpVar <- group_variable_name(moanin_model)
    tpVar <- time_variable_name(moanin_model)
    groups <- levels(droplevels(group_variable(moanin_model)) )
    
    # Check that the moanin model has the appropriate information to create a
    # smooth prediction model.
    if(is.null(spline_formula(moanin_model)) | 
       is.null(degrees_of_freedom(moanin_model))){
        msg <- paste(
            "Smooth prediction is not possible without the formula.",
            "Will only predict on initial points")
        warning(msg)
        return(moanin_model)
    }
    
    
    for(group in groups){
        mask <- group_variable(moanin_model) == group
        time <- time_variable(moanin_model)[mask]
        
        timepoints_pred <- c(
            timepoints_pred,
            seq(min(time), max(time), length=100))
        groups_pred <- c(
            groups_pred,
            rep(group, 100))
    }
    meta_prediction <- data.frame(
        "Timepoint"=timepoints_pred,
        "Replicates"=rep(1, length(timepoints_pred)),
        "Group"=groups_pred)
    names(meta_prediction)[c(1,3)]<-c(tpVar,gpVar)
    return(meta_prediction)
}


#' Rescales rows of data to be between 0 and 1
#'
#' @param data The matrix to rescale by row. If NULL, and \code{object} is
#'   given, data will be taken as \code{assay(object)} Each row should
#'   correspond to a gene or a centroid, and columns to samples.
#' @param object a object of class Moanin, only needed if choose to rescale by
#'   grouping variable in the moanin object. If NULL, then data will be rescaled
#'   jointly across all observations.
#' @param use_group If true, then the data will be rescaled such that, for each
#'   row, all values associated to each group (defined by grouping variable of
#'   \code{object}) is between 0 and 1. For example, if column
#' @param ... arguments passed to the matrix or Moanin method.
#' @details If the user set \code{log_transform=TRUE} in the creation of the
#'   \code{Moanin} object, the data will be log transformed before rescaling
#' @return rescaled y, such that for each row, the values are comprised between
#'   0 and 1. Note that if \code{use_group=TRUE} and \code{object} is not NULL,
#'   the values associated to the columns of unique values of the grouping
#'   variable of \code{object} will be rescaled separately.
#' @export
#' @name rescale_values
#' @examples 
#' data(exampleData)
#' moanin <- create_moanin_model(data=testData, meta=testMeta)
#' # Can rescale data in Moanin object
#' allData <- rescale_values(moanin)
#' # Or provide different data and/or rescale within grouping variable
#' smallData <- rescale_values(moanin, data=testData[1:10,], use_group=TRUE)
#' @aliases rescale_values,Moanin-method
setMethod("rescale_values", "Moanin",
    function(object, data=NULL, use_group=FALSE){
    if(is.null(data)){
        data <- get_log_data(object)
    }
    if(use_group){
        if(is.data.frame(data) || inherits(data, "DataFrame")){
            data <- data.matrix(data)
        }
        factors_to_consider <- levels(group_variable(object))
        for(factor in factors_to_consider){
            mask <- group_variable(object) == factor
            ymin <- matrixStats::rowMins(data[, mask]) 
            data[, mask] <- data[, mask] - ymin
            ymax <- matrixStats::rowMaxs(data[, mask])
            whNonZero <- which(ymax>0)
            data[whNonZero,mask] <- data[whNonZero,mask] / ymax[whNonZero]
        }
        return(data)
    } else {
        return(rescale_values(object=NULL, data=data))
    }
})

#' @aliases rescale_values, NULL-method
#' @export
#' @rdname rescale_values
setMethod("rescale_values", "NULL",
    function(object, data){
        if(inherits(data, "DataFrame")){
            data <- data.matrix(data)
        }
        ymin <- matrixStats::rowMins(data) 
        data <- data - ymin
        ymax <- matrixStats::rowMaxs(data)
        whNonZero <- which(ymax > 0)
        if(length(whNonZero) > 0){
            data[whNonZero,] <- data[whNonZero,] / ymax[whNonZero]
        }
        return(data)
})

#' @aliases rescale_values,missing-method
#' @export
#' @rdname rescale_values
setMethod("rescale_values", "missing",
          function(object, ...){
              rescale_values(object=NULL, ... )
          })

# Note that if would estimate a negative scaling factor, sets it to zero 
# instead, so that your centroid estimate for the data is a flat line at the mean
align_data_onto_centroid <- function(data, centroid, positive_scaling=TRUE, 
            returnType=c("data","centroid")){
    returnType <- match.arg(returnType)
    n_samples <- dim(data)[2]
    n_genes <- dim(data)[1]
    if(n_samples != length(centroid)){
        stop(
            paste0("align_data_onto_centroid: problem in dimensions. There are ",
                   n_samples, "but centroid is of length", length(centroid)))
    }
    
    # No clue why sometimes the vector/matrix is not numeric
    if(!is.numeric(centroid)){
        centroid <- as.numeric(centroid)
    }
    centered_centroid <- centroid - mean(centroid)
    # Identify if some rows of the data are only zeros.
    only_zero_genes <- matrixStats::rowSums2(abs(data)) == 0
    scaling_factors <- apply(
        data, 1,
        function(x){sum(centered_centroid * x)/sum((x - mean(x))*x)}) 
    if(positive_scaling){
        if(returnType=="data"){
            scaling_factors[scaling_factors < 0] <- 0
        }else{
            scaling_factors <- abs(scaling_factors)
        }
    }
    
    # Now replace the scaling factors of only 0 genes by 0
    scaling_factors[only_zero_genes] <- 0
    
    # Estimate the shift factors now
    shift_factors <- apply(
        scaling_factors * data,
        1,
        function(x) mean(centroid - x)) 
    if(returnType == "data"){
        data_fitted <- (scaling_factors * data + shift_factors)        
    }
    else{
        data_fitted <- matrix(centroid,
                              nrow=nrow(data),
                              ncol=length(centroid),
                              byrow=TRUE)
        data_fitted <- sweep(data_fitted,1,shift_factors,"-")
        whZero <- which(scaling_factors == 0)
        if(length(whZero) > 0){
            # If scaling factor is zero, then centroid estimate is just a flat 
            # line equal to the mean of the centroid - shift factor. 
            data_fitted[-whZero,] <- sweep(data_fitted[-whZero,,drop=FALSE],1,
                                          scaling_factors[-whZero],"/")
            data_fitted[whZero,] <- matrix(
                matrixStats::rowMeans2(data[whZero,,drop=FALSE]),
                nrow=length(whZero),
                ncol=ncol(data_fitted), byrow=FALSE)         
        }
        else 
            data_fitted <- sweep(data_fitted, 1, scaling_factors, "/")
        
    }
    return(data_fitted)
}


score_genes_centroid <- function(data, centroid, positive_scaling=TRUE, 
                                 scale=TRUE){
    n_genes <- dim(data)[1]
    centroid <- as.numeric(centroid)
    
    data_fitted <- align_data_onto_centroid(
        data, centroid, positive_scaling=positive_scaling)
    
    scores <- apply(
        data_fitted,
        1,
        function(y) sum((centroid - y)**2))
    
    if(scale){
        all_zeros_gene <- matrix(0, 1, dim(data)[2])
        all_zeros_gene <- align_data_onto_centroid(
            all_zeros_gene,
            centroid,
            positive_scaling=positive_scaling)
        score <- sum((centroid - all_zeros_gene)**2)
        max_score <- score
    }else{
        max_score <- 1
    }
    return(scores / max_score)
}


#' Fisher's method to combine pvalues
#'
#' Combines all p-values per rows.
#' 
#' @param pvalues a matrix of pvalues, with columns corresponding to different
#'   tests or sources of p-values, and rows corresponding to the genes from
#'   which the p-values come.
#' @return a vector of p-values, one for each row of \code{pvalues}, that is the
#'   result of Fisher's combined probability test applied to the p-values in
#'   that row.
#' @examples 
#' data(exampleData)
#' moanin <- create_moanin_model(data=testData,meta=testMeta)
#' contrasts <- create_timepoints_contrasts(moanin,"C", "K")
#' deTimepoints=DE_timepoints(moanin, 
#'   contrasts=contrasts, use_voom_weights=FALSE)
#' fisherPval=pvalues_fisher_method(
#'   deTimepoints[,grep("pval",colnames(deTimepoints))])
#' head(fisherPval)
#' @export
pvalues_fisher_method <- function(pvalues){
    # TODO Add a check that all pvalues are "valid"
    keep <- (pvalues >= 0) & (pvalues <= 1)
    pvalues[pvalues == 0] <- 1e-285
    
    lnp <- data.matrix(log(pvalues))
    chisq <- (-2) * matrixStats::rowSums2(lnp)
    df <- 2 * dim(lnp)[2]
    fisher_pval <- stats::pchisq(chisq, df, lower.tail=FALSE)
    return(fisher_pval)
}
