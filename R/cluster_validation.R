library(reshape2)
library(zoo)
library(NMI)


#' Compute consensus matrix from labels
#'
#' @param labels matrix (n_genes, n_clusters)
#'	Matrix containing the set of labels
#' @param scale	boolean, optional, default: TRUE
#'	Whether to rescale the matrix
#' @export
consensus_matrix = function(labels, scale=TRUE){
    melted_labels = reshape2::melt(as.matrix(labels))
    colnames(melted_labels) = c("Gene", "Clustering", "Label")
    melted_labels$Clustering_lab = melted_labels$Clustering:as.factor(melted_labels$Label)
    w = reshape2::dcast(melted_labels, Gene~Clustering_lab, value.var="Clustering_lab")
    x = as.matrix(w[,-1])
    x[is.na(x)] = 0
    x = apply(x, 2,  function(x) as.numeric(x > 0))
    consensus = tcrossprod(x) 
    if(scale){
        diag_elements = diag(consensus)
	consensus = t(t(consensus / diag_elements) / diag_elements)
	diag(consensus) = 0
    }
    return(consensus)
}

#' Plot CDF consensus matrix
#'
#' @param labels list of labels
#'
#' @export
plot_cdf_consensus = function(labels){
    all_labels = labels
    n_clusters = names(all_labels)

    colors = grDevices::rainbow(length(n_clusters))
    xrange = c(0, 1)
    yrange = c(0, 1)

    graphics::plot(xrange, yrange, type="n")

    for(i in 1:length(n_clusters)){
	cluster = n_clusters[i]
	color = colors[i]

	labels = all_labels[[cluster]]
	consensus = consensus_matrix(labels, scale=FALSE)
	consensus = consensus / max(consensus)
	consensus = sort(consensus[upper.tri(consensus)])
	x_axis = 1:length(consensus) / length(consensus)

	graphics::lines(consensus, x_axis, type="b",
			pch=16,
			col=color,
			lwd=1)

    }

    graphics::legend("bottomright", legend=n_clusters,
	lty=rep(1, length(n_clusters)),
       pch=rep(16, length(n_clusters)),
       col=colors,
       title="Clusters", text.font=4)

}

#' Plot AUC consensus matrix
#'
#' @param labels list of lables
#'
#' @export
get_auc_similarity_scores = function(labels, method="consensus"){
    if(!(method %in% c("consensus", "nmi"))){
	msg = paste(
	    "moanin::get_auc_similarity_score.",
	    "Unknown similarity measure provided", method)
	stop(msg)
    }

    all_labels = labels
    n_clusters = names(all_labels)

    auc_scores = rep(0, length(n_clusters))
    for(i in 1:length(n_clusters)){
	cluster = n_clusters[i]
	labels = all_labels[[cluster]]
	if(method == "consensus"){
	    consensus = consensus_matrix(labels, scale=FALSE)
	    consensus = consensus / max(consensus)
	    scores = consensus[upper.tri(consensus)]
	}else if(method == "nmi"){
	    scores = get_nmi_scores(labels) 
	}
	scores = sort(scores)
	y_axis = 1:length(scores) / length(scores)
	auc_score = sum(diff(scores) * zoo::rollmean(y_axis, 2))
	
	auc_scores[i] = auc_score
    }
    return(auc_scores)
}


get_nmi_scores = function(labels){

    nmi = function(x, y){
	x_dataframe = as.data.frame(x)
	colnames(x_dataframe) = c("Label")
	x_dataframe$Gene = row.names(x_dataframe)
	x_dataframe = x_dataframe[c("Gene", "Label")]

	y_dataframe = as.data.frame(y)
	colnames(y_dataframe) = c("Label")
	y_dataframe$Gene = row.names(y_dataframe)
	y_dataframe = y_dataframe[c("Gene", "Label")]

	return(NMI::NMI(x_dataframe, y_dataframe))
    }

    n_trials = dim(labels)[2]
    scores = NULL
    for(trial in 1:n_trials){
	if(trial == n_trials){
	    break
	}
	column = colnames(labels)[trial]
	columns_to_consider = colnames(labels)[(trial+1):n_trials]
	label = labels[column]
	scores = c(scores, as.vector(
	    unlist(apply(labels[columns_to_consider], 2, function(x){nmi(x, label)}))))
    }
    return(scores)
}


#' Plot model explorer
#'
#' @param labels list of labels
#' @export
plot_model_explorer = function(labels){
    all_labels = labels
    n_clusters = names(all_labels)

    nmi_scores = list()
    colors = grDevices::rainbow(length(n_clusters))
    max_trial = 0
    min_score = 1
    max_score = 0
    for(i in 1:length(n_clusters)){
	n_cluster = n_clusters[i]
	color = colors[i]

	labels = all_labels[[n_cluster]]
	scores = get_nmi_scores(labels)

	nmi_scores[[n_cluster]] = sort(scores)
	max_trial = max(max_trial, length(scores))
	min_score = min(min_score, min(scores))
	max_score = max(max_score, max(scores))
     }

    xrange = c(min_score, max_score)
    yrange = c(1, max_trial)

    graphics::plot(xrange, yrange, type="n", xlab="NMI", ylab="")

     for(i in 1:length(n_clusters)){
	color = colors[i]
	scores = nmi_scores[[i]]
	graphics::lines(sort(scores), 1:length(scores), type="b",
			pch=16,
			col=color,
			lwd=1)

    }

    graphics::legend(min_score, length(scores)*0.9, legend=n_clusters,
       lty=rep(1, length(scores)),
       pch=rep(16, length(scores)),
       col=colors, 
       title="Clusters", text.font=4)
}

