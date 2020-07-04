
#' Fit splines to each gene of data matrix
#'
#' @inheritParams DE_timecourse
#' @param weights A matrix of weights, of the same dimension as \code{data}.
#'
#' @return matrix of the coefficients for each basis function, each row of the
#'   matrix containing the coefficients for the corresponding gene in
#'   \code{data}.
#'
#' @export
fit_splines = function(data, moanin_model, weights=NULL){
    basis = moanin_model$basis
    n = ncol(basis)
    nr = nrow(data)
    
    if(!is.null(weights)){
        beta = matrix(nrow=nr, ncol=n)
        for(i in 1:nr){
            beta[i,] = stats::lm.wfit(basis, data[i,], weights[i,])$coefficients
        }
        row.names(beta) = row.names(data)
    }else{
        beta = t(stats::lm.fit(basis, t(data))$coefficients)
    }
    return(beta)
}

#' Get fitted values for splines
#'
#' @inheritParams DE_timecourse
#' @param meta_prediction optional, see \code{\link{create_meta_prediction}}.
#'
#' @return a matrix of the fitted y values, with dimensions the same as \code{data}
#'
#' @keywords internal
fit_predict_splines = function(data, moanin_model, 
                               meta_prediction=NULL){
    basis = moanin_model$basis
    meta = moanin_model$meta
    gpVar = moanin_model$group_variable
    tpVar = moanin_model$time_variable
    # if(!is.null(weights)){
    #     stop("moanin::fit_predict_splines: not implemented")
    # }
    if(is.null(meta_prediction)){
        y_fitted = t(stats::lm.fit(basis, t(data))$fitted.values)
    }else{
        degrees_of_freedom = moanin_model$degrees_of_freedom
        fitting_data = t(as.matrix(data))
        formula_data = list(
            "Group"=meta[,gpVar],
            "Timepoint"=meta[,tpVar],
            "fitting_data"=fitting_data,
            "degrees_of_freedom"=moanin_model$degrees_of_freedom)
        
        updated_formula = stats::update(moanin_model$formula, fitting_data ~ .)
        model = stats::lm(updated_formula, formula_data)	
        y_fitted = stats::predict(model, meta_prediction)
    }
    return(y_fitted)
}


#' Create prediction meta data from splines model
#'
#' @inheritParams DE_timecourse
#' @param num_timepoints integer, optional, default: 100.
#'	Number of timepoints to use for the prediction metadata
#'
#' @keywords internal
#' @export
create_meta_prediction = function(moanin_model, num_timepoints=100){
    # Create moanin_model for prediction
    timepoints_pred = NULL
    groups_pred = NULL
    meta = droplevels(moanin_model$meta)
    gpVar = moanin_model$group_variable
    tpVar = moanin_model$time_variable
    groups = levels(meta[,gpVar]) 

    # Check that the moanin model has the appropriate information to create a
    # smooth prediction model.
    if(is.null(moanin_model$formula) | is.null(moanin_model$degrees_of_freedom)){
        msg = paste(
            "Smooth prediction is not possible without the formula.",
            "Will only predict on initial points")
        warning(msg)
        return(moanin_model$meta)
    }


    for(group in groups){
    	mask = meta[,gpVar] == group
	    time = meta[,tpVar][mask]

    	timepoints_pred = c(
	        timepoints_pred,
	        seq(min(time), max(time), length=100))
	    groups_pred = c(
	        groups_pred,
	        rep(group, 100))
    }
    meta_prediction = data.frame(
	    "Timepoint"=timepoints_pred,
	    "Replicates"=rep(1, length(timepoints_pred)),
	    "Group"=groups_pred)
    names(meta_prediction)[c(1,3)]<-c(gpVar,tpVar)
    return(meta_prediction)
}


#' Rescales centroids and gene expresion values
#'
#' @param y 
#'      The matrix to rescale. Each row should correspond to a gene or a
#'      centroid and columns to samples.
#' @param meta, optional
#'      Metadata data.frame.
#' @param group, optional, default: NULL
#'      A column name of the metadata data.frame. The corresponding column
#'      should be factors. If provided, the values of y will be rescaled such
#'      that, for each row, all values associated to group A, … of column
#'      "group" of the metadata is between 0 and 1. For example, if column
#'      "group" corresponds to a genotype, all the values of a gene for a
#'      specific genotype will be rescaled between 0 and 1.
#' @return rescaled y, such that for each row, the values are comprised
#'      between 0 and 1. Note that if "group" is provided, the values
#'      associated to the columns of unique values of "group" will be rescaled
#'      separately.
#' @export
rescale_values = function(y, meta=NULL, group=NULL){
    if(is.null(group)){
	    ymin = row_min(y) 
	    y = y - ymin
	    ymax = row_max(y)
	    # We may have a division by 0 here
	    y = y / ymax
    }else{
        if(is.null(meta)){
            msg = paste(
                "moanin::rescale_values if group is provided, then a metadata",
                "data.frame should be provided as well.")
            stop(msg)
        }
    	factors_to_consider = levels(unlist(meta[group]))
	    for(factor in factors_to_consider){
	        mask = meta[group] == factor
            ymin = row_min(y[, mask]) 
            y[, mask] = y[, mask] - ymin
            ymax = row_max(y[, mask])
            # We may have a division by 0 here
            y[, mask] = y[, mask] / ymax
        }
    }
    return(y)
}


# XXX It's wierd that this does not exists in R…
# it probably exists but under another name?
row_max = function(X){
  return(apply(X, 1, max))
}

row_min = function(X){
  return(apply(X, 1, min))
}

row_mean = function(X){
    return(apply(X, 1, mean))
}

row_sum = function(X){
    return(apply(X, 1, sum))
}

row_argmin = function(X){
    return(apply(X, 1, which.min))
}


# Worst name ever
align_data_onto_centroid = function(data, centroid, positive_scaling=TRUE){
    n_samples = dim(data)[2]
    n_genes = dim(data)[1]
    if(n_samples != length(centroid)){
	stop("align_data_onto_centroid: problem in dimensions")
    }

    # No clue why sometimes the vector/matrix is not numeric
    if(!is.numeric(centroid)){
	centroid = as.numeric(centroid)
    }
    centered_centroid = centroid - mean(centroid)
    # Identify if some rows of the data are only zeros.
    only_zero_genes = rowSums(abs(data)) == 0
    scaling_factors = apply(
	data, 1,
	function(x){sum(centered_centroid * x)/sum((x - mean(x))*x)}) 
    if(positive_scaling){
        scaling_factors[scaling_factors < 0] = 0
    }

    # Now replace the scaling factors of only 0 genes by 0
    scaling_factors[only_zero_genes] = 0

    # Estimate the shift factors now
    shift_factors = apply(
	scaling_factors * data,
	1,
	function(x) mean(centroid - x)) 

    data_fitted = (scaling_factors * data + shift_factors)
    return(data_fitted)
}


score_genes_centroid = function(data, centroid, positive_scaling=TRUE, scale=TRUE){
    n_genes = dim(data)[1]
    centroid = as.numeric(centroid)

    data_fitted = align_data_onto_centroid(
	data, centroid, positive_scaling=positive_scaling)

    scores = apply(
	data_fitted,
	1,
	function(y) sum((centroid - y)**2))

    if(scale){
	all_zeros_gene = matrix(0, 1, dim(data)[2])
        all_zeros_gene = align_data_onto_centroid(
	    all_zeros_gene,
	    centroid,
	    positive_scaling=positive_scaling)
	score = sum((centroid - all_zeros_gene)**2)
	max_score = score
    }else{
	max_score = 1
    }
    return(scores / max_score)
}


#' Fisher's method to combine pvalues
#'
#' Combines all p-value per rows.
#' 
#' @param pvalues a matrix of pvalues, with columns corresponding to different
#'   tests or sources of p-values, and rows corresponding to the genes from
#'   which the p-values come.
#' @return a vector of p-values, one for each row of \code{pvalues}, that is the
#'   result of Fisher's combined probability test applied to the p-values in
#'   that row.
#' @export
pvalues_fisher_method = function(pvalues){
    # TODO Add a check that all pvalues are "valid"
    keep = (pvalues >= 0) & (pvalues <= 1)
    pvalues[pvalues == 0] = 1e-285
    
    lnp = log(pvalues)
    chisq = (-2) * row_sum(lnp)
    df = 2 * length(lnp)
    fisher_pval = stats::pchisq(chisq, df, lower.tail=FALSE)
    return(fisher_pval)
}
