library(ClusterR)
library(stats)
library(splines)

#' Performs splines clustering using K-means
#'
#' @param data matrix containg the data. Data needs to be in log scale.
#' @param moanin_model moanin_model 
#' @param n_clusters int optional, default: 10
#' @param init	["kmeans++", "random", "optimal_init"]
#' @param n_init int, optional, default: 10
#'	    Number of initialization to perform.
#' @param max_iter  int, optional, default: 300
#'	Maximum number of iteration to perform	
#' @param random_seed int, optional, default: NULL
#' @param fit_splines	boolean, optional, default: TRUE
#'	Whether to fit splines or not.
#' @param rescale   boolean, optional, default: TRUE
#'	Whether to rescale the data or not.
#' @export
splines_kmeans = function(data, moanin_model, n_clusters=10,
			  init="kmeans++",
			  n_init=10,
			  max_iter=300,
			  random_seed=NULL,
			  fit_splines=TRUE,
			  rescale=TRUE){
    meta = moanin_model$meta
    basis = moanin_model$basis
    check_data_meta(data, meta)

    if(fit_splines){
        fitted_data = fit_predict_splines(data, moanin_model)
    }else{
	fitted_data = data
    }

    if(rescale){
        fitted_data = rescale_values(fitted_data, meta)
    }

    # Set the random seed if it is null.
    if(is.null(random_seed)){
	set.seed(NULL)
	random_seed = .Random.seed[1]
    }
    kmeans_clusters = ClusterR::KMeans_rcpp(
	fitted_data, n_clusters, num_init=n_init, max_iters=max_iter,
	seed=random_seed, initializer=init)
    kmeans_clusters$centroids = rescale_values(
	kmeans_clusters$centroids, moanin_model)
    names(kmeans_clusters$clusters) = row.names(data)

    # Give names to clusters
    cluster_names = sapply(1:n_clusters, function(x){paste0("C", x)})
    row.names(kmeans_clusters$centroids) = cluster_names
    colnames(kmeans_clusters$centroids) = colnames(data)

    kmeans_clusters$moanin_model = moanin_model
    kmeans_clusters$fit_splines = fit_splines
    kmeans_clusters$rescale = rescale
    return(kmeans_clusters)
}


splines_kmeans_prediction = function(data, kmeans_clusters){
    moanin_model = kmeans_clusters$moanin_model
    fit_splines = kmeans_clusters$fit_splines
    rescale = kmeans_clusters$rescale

    meta = moanin_model$meta
    basis = moanin_model$basis
    check_data_meta(data, meta)

    if(fit_splines){
        fitted_data = fit_predict_splines(data, moanin_model)
    }else{
	fitted_data = data
    }

    if(rescale){
        fitted_data = rescale_values(fitted_data, meta)
    }

    closest_cluster <- function(x) {
	cluster_dist <- apply(
	    kmeans_clusters$centroids, 1, function(y){sqrt(sum((x-y)^2))})
	return(which.min(cluster_dist)[1])
    }

    all_labels <- apply(fitted_data, 1, closest_cluster) 
    kmeans_clusters$clusters = all_labels
    names(kmeans_clusters$clusters) = row.names(data)
    return(kmeans_clusters)
}

#' Assign score and labels from raw data
#'
#' @param data array or dataframe.
#'	Data to assign score and labels to
#' @param kmeans_clusters list of list
#'	List returned by moanin::splines_kmeans
#' @param percentage_genes_to_label float, optional, default: 0.5
#'	Percentage of genes to label. If max_score is provided, will label
#'	genes that are either in the top `percentage_genes_to_label` or with a
#'	score below `max_score`.
#' @param max_score optional, default: Null
#'	When provided, will only label genes below that score. If NULL, ignore
#'	this option.
#' @param rescale_separately_on, string, optional, default: NULL
#'	When provided, will rescale separately different groups of data.
#' @export
splines_kmeans_score_and_label = function(data, kmeans_clusters, percentage_genes_to_label=0.5,
					  max_score=NULL,
					  rescale_separately_on=NULL){

    n_clusters = dim(kmeans_clusters$centroids)[1]
    all_scores = matrix(NA, nrow=dim(data)[1], ncol=n_clusters)


    if(!is.null(rescale_separately_on)){
        groups = levels(meta[, rescale_separately_on])
    }

    for(k in 1:n_clusters){
	# By default, should not rescale separately on any columns.
	if(is.null(rescale_separately_on)){
	
	    scores = score_genes_centroid(
		data,
		kmeans_clusters$centroids[k,],
		scale=FALSE)

	    all_scores[, k] = scores /  max(scores)
	}else{
	    scores = NULL
	    for(group in groups){
		mask = meta[, rescale_separately_on] == group
		partial_scores = score_genes_centroid(
		    data[, mask],
		    kmeans_clusters$centroids[k, mask],
		    scale=FALSE)
		if(is.null(scores)){
		    scores = partial_scores
		}else{
		    scores = scores + partial_scores
		}
	    }	
	    all_scores[, k] = scores / max(scores)
	}
    }

    # Give names to rows
    all_scores = as.matrix(all_scores)
    row.names(all_scores) = row.names(data) 
    scores = apply(all_scores, 1, min)

    # Only assign labels to X% of the genes
    max_score_data = stats::quantile(scores, c(percentage_genes_to_label))
    if(!is.null(max_score)){
	max_score = min(max_score_data, max_score)
    }else{
	max_score = max_score_data
    }
    genes_to_not_consider = scores >= max_score
    labels = apply(all_scores, 1, which.min)
    labels[genes_to_not_consider] = NA
    names(labels) = row.names(data)

    if(max_score == 1){
	msg = paste(
	    "moanin::splines_kmeans_score_and_label is labeling genes",
	    " with a score of 1. This implies that no good fit for ",
	    "those genes are found and the assignement is random.",
	    sep="")
	warning(msg)
    }

    return(list("labels"=labels, "scores"=all_scores))
}
