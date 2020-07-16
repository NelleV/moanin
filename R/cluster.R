
setGeneric("splines_kmeans", 
           function(object,...) { standardGeneric("splines_kmeans")})
setGeneric("splines_kmeans_score_and_label", 
           function(object,...) { standardGeneric("splines_kmeans_score_and_label")})
setGeneric("splines_kmeans_predict", 
           function(object,...) { standardGeneric("splines_kmeans_predict")})


#' Performs splines clustering using K-means
#'
#' @inheritParams DE_timecourse
#' @param n_clusters int optional, default: 10
#' @param init  ["kmeans++", "random", "optimal_init"]
#' @param n_init int, optional, default: 10
#'      Number of initialization to perform.
#' @param max_iter  int, optional, default: 300
#'  Maximum number of iteration to perform  
#' @param random_seed int, optional, default: NULL. 
#' Passed to argument \code{seed} in  \code{\link[ClusterR]{KMeans_rcpp}}. 
#' If NULL (default), set to \code{.Random.seed[1]}.
#' @param fit_splines   boolean, optional, default: TRUE
#'  Whether to fit splines or not.
#' @param rescale   boolean, optional, default: TRUE
#'  Whether to rescale the data or not.
#' @details If \code{\link{Moanin}} object's slot has log_transform=TRUE, then the data will be transformed by the function log(x+1) before applying splines and clustering.
#' @return A list in the format returned by \code{\link[ClusterR]{KMeans_rcpp}},
#'    with the following elements added or changed:
#'\itemize{
#'\item{\code{centroids}}{ The centroids are rescaled so that they range from
#'0-1}
#'\item{\code{fit_splines}}{ Logical, the value of \code{fit_splines} given to the
#'function }
#'\item{\code{rescale}}{ The value of \code{rescale} given to the function }
#'} 
#' @name splines_kmeans
#' @aliases splines_kmeans,Moanin-method
#' @examples 
#' data(exampleData)
#' # Use the default options
#' moanin = create_moanin_model(data=testData, meta=testMeta)
#' out = splines_kmeans( moanin,n_clusters=5)
#' table(out$clusters)
#' @importFrom ClusterR KMeans_rcpp
#' @export
setMethod("splines_kmeans", "Moanin",
    function(object, n_clusters=10,
                        init="kmeans++",
                        n_init=10, 
                        max_iter=300, 
                        random_seed=.Random.seed[1],
                        fit_splines=TRUE,
                        rescale=TRUE){
    basis = basis_matrix(object)
    data<-get_log_data(object)
    if(fit_splines){
        fitted_data = fit_predict_splines(data=data,object)
    }else{
        fitted_data = data
    }
    
    ## CHECK ME: previous version gave the meta information to rescale_values, 
    ## but also group=NULL, so meant it wasn't used. Was that a mistake?
    if(rescale){
        fitted_data = rescale_values(data=fitted_data, object=object,
                                     use_group=FALSE)
    }
    
    kmeans_clusters = ClusterR::KMeans_rcpp(
        fitted_data, n_clusters, num_init=n_init, max_iters=max_iter,
        seed=random_seed, initializer=init)
    kmeans_clusters$centroids = rescale_values(
        data=kmeans_clusters$centroids, object=object,use_group=FALSE )
    names(kmeans_clusters$clusters) = row.names(object)
    
    # Give names to clusters
    cluster_names = vapply(seq_len(n_clusters), FUN=function(x){paste0("C", x)},
                        FUN.VALUE="C")
    row.names(kmeans_clusters$centroids) = cluster_names
    colnames(kmeans_clusters$centroids) = colnames(object)
    
    #kmeans_clusters$moanin_model = moanin_model
    kmeans_clusters$fit_splines = fit_splines
    kmeans_clusters$rescale = rescale
    return(kmeans_clusters)
}
)

#' @param method If "distance", predicts based on distance of data to kmeans
#'   centroids. If "goodnessOfFit", is a wrapper to
#'   \code{splines_kmeans_score_and_label}, assigning labels based on goodness
#'   of fit, including any filtering.
#' @param ... arguments passed to \code{splines_kmeans_score_and_label}
#' @return \code{splines_kmeans_predict} returns a vector giving the labels for
#'   the given data.
#' @rdname splines_kmeans_score_and_label
#' @export
setMethod("splines_kmeans_predict", "Moanin",
          function(object, kmeans_clusters,data=NULL, 
                   method=c("distance","goodnessOfFit"),...){
    method=match.arg(method)
    if(method=="goodnessOfFit"){
        out<-splines_kmeans_score_and_label(object=object, 
            kmeans_clusters=kmeans_clusters,data=data,...)
        return(out$labels)
    }
    if(method=="distance"){
        fit_splines = kmeans_clusters$fit_splines
        rescale = kmeans_clusters$rescale
        check_data_meta(kmeans_clusters$centroids, object)
        basis = basis_matrix(object)
        if(is.null(data)){data<-get_log_data(object)
        }
        else{
            check_data_meta(data,object)
        }
        if(fit_splines){
            fitted_data = fit_predict_splines(data=data, moanin_model=object)
        }else{
            fitted_data = data
        }
        
        if(rescale){
            fitted_data = rescale_values(data=fitted_data, object=object, 
                                         use_group=FALSE)
        }
        
        closest_cluster <- function(x) {
            cluster_dist <- apply(
                kmeans_clusters$centroids, 1, function(y){sqrt(sum((x-y)^2))})
            return(which.min(cluster_dist)[1])
        }
        
        all_labels <- apply(fitted_data, 1, closest_cluster) 
        kmeans_clusters$clusters = all_labels
        names(kmeans_clusters$clusters) = row.names(data)
        return(kmeans_clusters$clusters)
    }

}
)


#' Assign score and labels from raw data
#' @name splines_kmeans_score_and_label
#' @param object the Moanin object that contains the basis functions used in
#'   creating the clusters
#' @param kmeans_clusters the results of running \code{\link{splines_kmeans}}
#' @param data the data to predict. If not given, will use \code{assay(object)}.
#'   If given, the number of columns of \code{data} must match that of
#'   \code{object}
#' @param kmeans_clusters  List returned by \code{\link{splines_kmeans}}
#' @param proportion_genes_to_label float, optional, default: 0.5
#'  Percentage of genes to label. If max_score is provided, will label
#'  genes that are either in the top `proportion_genes_to_label` or with a
#'  score below `max_score`.
#' @param max_score optional, default: Null
#'  When provided, will only label genes below that score. If NULL, ignore
#'  this option.
#' @param previous_scores matrix of scores, optional. Allows user to give the
#'   matrix scores results from a previous run of
#'   \code{splines_kmeans_score_and_label}, and only redo the filtering (i.e. if
#'   want to change \code{proportion_genes_to_label} without rerunning the
#'   calculation of scores)
#' @param rescale_separately logical, whether to score separately within
#'   grouping variable
#' @return A list consisting of
#' \itemize{
#' \item{\code{labels}}{the label or cluster assigned to each gene based on the
#' cluster with the best (i.e. lowest) score, with no label given to genes that
#' do not have a score lower than a specified quantity}
#' \item{\code{scores}}{the matrix of size n_cluster x n_genes, containing for 
#' each gene and each cluster, the goodness of fit score}
#' \item{\code{score_cutoff}}{The required cutoff for a gene receiving an
#' assignment}
#' }
#' @aliases splines_kmeans_predict splines_kmeans_predict,Moanin-method 
#' @aliases splines_kmeans_score_and_label,Moanin-method
#' @examples 
#' data(exampleData)
#' moanin = create_moanin_model(data=testData, meta=testMeta)
#' # Cluster on a subset of genes
#' kmClusters=splines_kmeans(moanin[1:50,],n_clusters=3)
#' # get scores on all genes
#' scores_and_labels = splines_kmeans_score_and_label(object=moanin, kmClusters)
#' head(scores_and_labels$scores)
#' head(scores_and_labels$labels)
#' # should be same as above, only just the assignments
#' predictLabels1 = splines_kmeans_predict(object=moanin, kmClusters, 
#'      method="goodnessOfFit")
#' # Instead use distance to centroid:
#' predictLabels2 = splines_kmeans_predict(object=moanin, kmClusters, 
#'      method="distance")
#' @export
setMethod("splines_kmeans_score_and_label", "Moanin",
          function(object, kmeans_clusters, data=NULL,
                   proportion_genes_to_label=0.5,
                   max_score=NULL, previous_scores=NULL,
                   rescale_separately=FALSE){
    if(is.null(data)){
        data<-get_log_data(object)
    }
    else{
        if(ncol(data)!=ncol(kmeans_clusters$centroids)) stop(
            "User-given data and kmeans resultsare inconsistent. Data is has ", 
            ncol(data),"columns; kmeans result was run on", 
            ncol(kmeans_clusters$centroids), "samples")
    }

    if(is.null(previous_scores)){
        n_clusters = dim(kmeans_clusters$centroids)[1]
        all_scores = matrix(NA, nrow=dim(data)[1], ncol=n_clusters)


        for(k in seq_len(n_clusters)){
            # By default, should not rescale separately on any columns.
            if(!rescale_separately){
                scores = score_genes_centroid(
                    data,
                    kmeans_clusters$centroids[k,],
                    scale=FALSE)

                all_scores[, k] = scores /  max(scores)
            }else{
                scores = NULL
                groups = levels(group_variable(object))
                for(group in groups){
                    mask = group_variable(object) == group
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
    }
    else all_scores=previous_scores
    
    scores = apply(all_scores, 1, min)
    labels = apply(all_scores, 1, which.min)
    names(labels) = row.names(data)

    if(proportion_genes_to_label<1 | !is.null(max_score)){
        max_score_data = stats::quantile(scores, c(proportion_genes_to_label))
        if(!is.null(max_score)){
            max_score = min(max_score_data, max_score)
        }else{
            max_score = max_score_data
        }
        genes_to_not_consider = scores >= max_score
        labels[genes_to_not_consider] = NA
    
        if(max_score == 1){
            msg = paste(
                "moanin::splines_kmeans_score_and_label is labeling genes",
                " with a score of 1. This implies that no good fit for ",
                "those genes are found and the assignment is random.",
                sep="")
            warning(msg)
        } 
    }# Only assign labels to X% of the genes


    return(list("labels"=labels,
                "scores"=all_scores,
                "score_cutoff"=max_score))
          }
)
