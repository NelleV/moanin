library(viridis)
library(graphics)


#' Plotting centroids
#'
#' @param centroids matrix (k, t) containing the centroids
#' @param splines_model splines_model
#' @param meta	data.frame (t, n) containing the metadata.
#' @param colors vector, optional, default NULL
#'		vector of colors
#' @param smooth boolean, optional, default: FALSE
#'  Whether to smooth the centroids or not.
#' @param margins vector containing the margins between
#' @export
plot_centroids = function(centroids, splines_model, colors=NULL, smooth=FALSE, mar=c(1.0, 1.0, 1.0, 1.0)){
    n_centroids = dim(centroids)[1]
    if(n_centroids <= 3){
        graphics::par(
	    mfrow=c(n_centroids, 1),
            mar=c(1.0, 1.0, 1.0, 1.0))
    }else if(n_centroids <= 6){
        n_col = ceiling(n_centroids / 2)
        graphics::par(
	    mfrow=c(n_col, 2),
            mar=c(1.0, 1.0, 1.0, 1.0))
    }else if(n_centroids <= 12){
        ncol = ceiling(n_centroids / 3)
        graphics::par(
	    mfrow=c(ncol, 3),
            mar=c(1.0, 1.0, 1.0, 1.0))
    }else{
        nrow = round(n_centroids ** 0.5)
        ncol = ceiling(n_centroids / nrow)
        graphics::par(
	    mfrow=c(ncol, nrow),
            mar=c(1.0, 1.0, 1.0, 1.0))
    }

    name_centroids = row.names(centroids)
    name_centroid = NULL
    for(i in 1:n_centroids){
	if(!is.null(name_centroids)){
	    name_centroid = name_centroids[i]
	}
        plot_centroid_individual(as.vector(centroids[i, ]),
				 splines_model, colors=colors,
				 smooth=smooth,
				 title=name_centroid)
    }
}


#' Plotting data
#'
#' @param data matrix (k, t) containing the centroids
#' @param splines_model splines_model
#' @param meta	data.frame (t, n) containing the metadata.
#' @param colors vector, optional, default NULL
#'		vector of colors
#' @param smooth boolean, optional, default: FALSE
#'  Whether to smooth the centroids or not.
#' @export
plot_genes = function(data, splines_model, colors=NULL, smooth=FALSE){
    plot_centroids(data, splines_model, colors=colors, smooth=smooth)
}


plot_centroid_individual = function(centroid, splines_model, colors=NULL, smooth=FALSE,
				    title=NULL){
    meta = splines_model$meta
    groups = levels(meta$Group)

    xrange = range(meta$Timepoint)
    yrange = range(centroid)
    if(is.null(dim(centroid))){
        centroid = t(as.matrix(centroid))
    }

    graphics::plot(xrange, yrange, type="n", main=title)
    if(is.null(colors)){
        colors = viridis::viridis(length(groups))
	names(colors) = groups
    }
    if(smooth){
	# FIXME this is supposed to be on the fitted lines, but I'm not able
	# to get this to work fine in R.
	meta_prediction = create_meta_prediction(splines_model)
	centroid_fitted = fit_predict_splines(
	    centroid, splines_model,
	    meta_prediction=meta_prediction)
    }else{
	centroid_fitted = fit_predict_splines(
	    centroid, splines_model)
    }


    # scatter points for values
    for(i in 1:length(groups)){
        group = groups[i]
        color = colors[group]
        
	# Start by individual points
        mask = meta$Group == group
        time = meta$Timepoint[mask]
        indx = order(time)
        graphics::lines(time[indx], centroid[mask][indx], type="p",
		        col=color, pch=16,
			lwd=0)
	if(smooth){
	    mask = meta_prediction$Group == group
	    time = meta_prediction$Timepoint[mask]
	    indx = order(time)
	}
	graphics::lines(time[indx], centroid_fitted[mask][indx], type="l",
			col=color, lwd=1)

    }
} 

plot_gene_splines = function(data, meta, gene_name, colors=NULL){
    # First, select the gene:
    if(!(gene_name %in% row.names(data))){
	msg = paste("moanin::plot_gene_splines: The gene_name provided '",
		    gene_name, "' is not in the data", sep="")
    }
    gene_data = data[gene_name, ]


}
