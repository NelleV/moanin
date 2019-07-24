library(ClusterR)
library(stats)
library(splines)

#' Performs splines clustering using K-means
#'
#' @param data matrix containg the data. Data needs to be in log scale.
#' @param splines_model splines_model 
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
splines_kmeans = function(data, splines_model, n_clusters=10,
			  init="kmeans++",
			  n_init=10,
			  max_iter=300,
			  random_seed=NULL,
			  fit_splines=TRUE,
			  rescale=TRUE){
    meta = splines_model$meta
    basis = splines_model$basis
    check_data_meta(data, meta)

    if(fit_splines){
        fitted_data = fit_predict_splines(data, splines_model)
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
	kmeans_clusters$centroids, splines_model)
    names(kmeans_clusters$clusters) = row.names(data)

    # Give names to clusters
    cluster_names = sapply(1:n_clusters, function(x){paste0("C", x)})
    row.names(kmeans_clusters$centroids) = cluster_names
    colnames(kmeans_clusters$centroids) = colnames(data)

    kmeans_clusters$splines_model = splines_model
    kmeans_clusters$fit_splines = fit_splines
    kmeans_clusters$rescale = rescale
    return(kmeans_clusters)
}


splines_kmeans_prediction = function(data, kmeans_clusters){
    splines_model = kmeans_clusters$splines_model
    fit_splines = kmeans_clusters$fit_splines
    rescale = kmeans_clusters$rescale

    meta = splines_model$meta
    basis = splines_model$basis
    check_data_meta(data, meta)

    if(fit_splines){
        fitted_data = fit_predict_splines(data, splines_model)
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
#'	Percentage of genes to label.
#' @export
splines_kmeans_score_and_label = function(data, kmeans_clusters, percentage_genes_to_label=0.5){

    n_clusters = dim(kmeans_clusters$centroids)[1]
    all_scores = matrix(NA, nrow=dim(data)[1], ncol=n_clusters)

    for(k in 1:n_clusters){
	scores = score_genes_centroid(
	    data,
	    kmeans_clusters$centroids[k,])
	all_scores[, k] = scores
    }

    # Give names to rows
    all_scores = as.matrix(all_scores)
    row.names(all_scores) = row.names(data) 

    # Only assign labels to X% of the genes
    max_score = stats::quantile(rowMin(all_scores), c(percentage_genes_to_label))
    genes_to_not_consider = rowMin(all_scores) >= max_score
    labels = apply(all_scores, 1, which.min)
    labels[genes_to_not_consider] = NA
    names(labels) = row.names(data)

    return(list("labels"=labels, "scores"=all_scores))
}
