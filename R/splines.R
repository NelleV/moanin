# Splines utility

library(splines)
library(MASS)
library(stats)


#' Fit splines
#'
#' @param data the data
#' @param splines_model splines_model 
#' @param weights weigts
#'
#' @return beta coefficients
#'
#' @export
fit_splines = function(data, splines_model, weights=NULL){
    basis = splines_model$basis
    n = ncol(basis)
    nr = nrow(data)
    basis = splines_model$basis 

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

#' Fit and predict splines
#'
#' @param data the data
#' @param splines_model splines_model
#' @param weights weights
#'
#' @return y_fitted the fitted y values
#'
#' @export
fit_predict_splines = function(data, splines_model, weights=NULL, meta_prediction=NULL){
    basis = splines_model$basis
    meta = splines_model$meta
    if(!is.null(weights)){
	stop("moanin::fit_predict_splines: not implemented")
    }
    if(is.null(meta_prediction)){
        y_fitted = t(stats::lm.fit(basis, t(data))$fitted.values)
    }else{
	degrees_of_freedom = splines_model$degrees_of_freedom
	fitting_data = t(as.matrix(data))
	formula_data = list(
	    "Group"=meta$Group,
	    "Timepoint"=meta$Timepoint,
	    "fitting_data"=fitting_data,
	    "degrees_of_freedom"=splines_model$degrees_of_freedom)

	updated_formula = stats::update(splines_model$formula, fitting_data ~ .)
	model = stats::lm(updated_formula, formula_data)	
	y_fitted = stats::predict(model, meta_prediction)
    }
    return(y_fitted)
}


#' Create prediction meta data from splines model
#'
#' @param splines_model a Splines model object
#' @param num_timepoints integer, optional, default: 100
#'	Number of timepoints to use for the prediction metadata
#'
#' @export
create_meta_prediction = function(splines_model, num_timepoints=100){
    # Create splines_model for prediction
    timepoints_pred = NULL
    groups_pred = NULL
    meta = droplevels(splines_model$meta)
    groups = levels(meta$Group) 
    for(group in groups){
	mask = meta$Group == group
	time = meta$Timepoint[mask]

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
    return(meta_prediction)
}

rescale_values = function(y, meta=NULL, group=NULL){
    if(is.null(group)){
	ymin = row_min(y) 
	y = y - ymin
	ymax = row_max(y)
	# We may have a division by 0 here
	y = y / ymax
    }else{
	factors_to_consider = levels(unlist(meta[group]))
	for(factor in factors_to_consider){
	    mask = meta["Group"] == factor
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
    shift_factors = rowMeans(
	rep(centroid, each=n_genes) - rep(scaling_factors, times=n_samples) * data)

    data_fitted = rep(scaling_factors, times=n_samples) * data
    data_fitted = (data_fitted + rep(shift_factors, times=n_samples))
    return(data_fitted)
}


score_genes_centroid = function(data, centroid, positive_scaling=TRUE, scale=TRUE){
    n_genes = dim(data)[1]
    data_fitted = align_data_onto_centroid(
	data, centroid, positive_scaling=positive_scaling)

    scores = apply(
	data_fitted,
	1,
	function(y){sqrt(sum((centroid - y)^2))})

    if(scale){
	all_zeros_gene = matrix(0, 1, dim(data)[2])
        all_zeros_gene = align_data_onto_centroid(
	    all_zeros_gene,
	    centroid,
	    positive_scaling=positive_scaling)
	score = sqrt(sum((centroid - all_zeros_gene)**2))
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
#' @param pvalues pvalues
#'
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
