# library(viridis)
# library(graphics)


#' Plotting centroids
#'
#' @param centroids matrix (k, t) containing the centroids
#' @param splines_model splines_model
#' @param meta	data.frame (t, n) containing the metadata.
#' @param colors vector, optional, default NULL
#'		vector of colors
#' @param smooth boolean, optional, default: FALSE
#'  Whether to smooth the centroids or not.
#' @param mar vector of margins to set the space around each plot (see \code{\link{par}})
#' @param legend boolean whether to include a legend (default:TRUE)
#' @param legendArgs list of arguments to be passed to legend command (if \code{legend=TRUE})
#' @param ... arguments to be passed to the individual plot commands
#'  (Will be sent to all plot commands)
#' @param simpleY boolean, if true, will minimize the annotation of the y axis to 
#'  only label the axis in the exterior plots (the x-axis is always assumed to be the 
#'  same across all plots and will always be simplified)
#' @param mfrow a vector of integers of length 2 defining the grid of plots to be created (see \code{\link{par}}). If missing, the function will set a value.

#' @export
plot_centroids = function(centroids, splines_model, colors=NULL, smooth=FALSE, legend=TRUE, legendArgs=NULL,simpleY=TRUE, mar=c(2.5, 2.5, 3.0, 1.0),mfrow=NULL,...){
    n_centroids = dim(centroids)[1]
    n_plots<-if(legend) n_centroids+1 else n_centroids
    if(!is.null(mfrow)){
        if(mfrow[1]*mfrow[2] < n_plots) stop(sprintf("Invalid value for argument mfrow. Should result in grid for at least %s plots (including a plot for the legend, if legend=TRUE)", n_plots))
    }
    else{
        if(n_plots <= 3){ 
            mfrow<-c(n_plots, 1)
        }else if(n_plots <= 6){
            mfrow=c(ceiling(n_plots / 2), 2)
        }else if(n_plots <= 12){
            mfrow=c(ceiling(n_plots / 3), 3)
        }else{
            nrow = round(n_plots ** 0.5)
            ncol = ceiling(n_plots / nrow)
            mfrow=c(ncol, nrow)
        }        
    }
    graphics::par(mfrow=mfrow, mar=mar)
    bottomPlots<-seq(to=n_centroids,by=1,length=mfrow[2])
    sidePlots<-seq(from=1,to=n_centroids,by=mfrow[2])

    name_centroids = row.names(centroids)
    name_centroid = NULL
    
    if(is.null(colors)){
        groups = levels(meta$Group)
        colors = viridis::viridis(length(groups))
        names(colors) = groups
    }
    for(i in 1:n_centroids){
    	if(!is.null(name_centroids)){
    	    name_centroid = name_centroids[i]
    	}
        plot_centroid_individual(as.vector(centroids[i, ]),
				 splines_model, colors=colors,
				 smooth=smooth,
				 title=name_centroid,
                 xaxt=if(!i %in% bottomPlots) "n" else "s",
                 yaxt=if(!i %in% sidePlots & simpleY) "n" else "s", ...
                 )
        if(!i %in% bottomPlots) axis(1, labels=FALSE)
        if(!i %in% sidePlots & simpleY) axis(2, labels=FALSE)
        
    }
    if(legend){
        plot.new()
        plot.window(xlim=c(0,1), ylim=c(0,1), bty="n",xaxt="n",yaxt="n")
        do.call("legend",c(list(x="center",
            legend=names(colors),fill=colors,bty="n"),legendArgs))
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
#' @param ... arguments passed to plot_centroids
#' @export
plot_genes = function(data, splines_model, colors=NULL, smooth=FALSE,...){
    plot_centroids(data, splines_model, colors=colors, smooth=smooth,simpleY=FALSE,...)
}


plot_centroid_individual = function(centroid, splines_model, colors, smooth=FALSE,
				    title=NULL,...){
    meta = splines_model$meta
    groups = levels(meta$Group)

    xrange = range(meta$Timepoint)
    yrange = range(centroid)
    if(is.null(dim(centroid))){
        centroid = t(as.matrix(centroid))
    }

    graphics::plot(xrange, yrange, type="n", main=title,...)
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
