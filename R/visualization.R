# library(viridis)
# library(graphics)


#' Plotting splines
#'
#' @param data matrix (k, t) containing the data (such as centroids or data) 
#' @param moanin_model moanin_model
#' @param colors vector, optional, default NULL
#'		vector of colors
#' @param smooth boolean, optional, default: FALSE
#'  Whether to smooth the centroids or not.
#' @param mar vector of margins to set the space around each plot
#'   (see \code{\link{par}})
#' @param legend boolean whether to include a legend (default:TRUE)
#' @param legendArgs list of arguments to be passed to legend command
#'    (if \code{legend=TRUE})
#' @param simpleY boolean, if true, will minimize the annotation of the
#'   y axis to  only label the axis in the exterior plots (the x-axis
#'   is always assumed to be the same across all plots and will always
#'    be simplified)
#' @param subset_conditions list
#'   if provided, only plots the subset of conditions provided. Else, plots
#'   all conditions
#' @param subset_data list
#'   if provided, only plots  the subset of data (ie, the rows) provided.
#' @param mfrow a vector of integers of length 2 defining the grid of
#'   plots to be created (see \code{\link{par}}). If missing, the
#'   function will set a value.
#' @param addToPlot A function that will be called after the plotting, allowing the user to add more to the plot.
#' @param ... arguments to be passed to the individual plot commands
#'  (Will be sent to all plot commands)
#'
#' @examples
#' # First, load some data and create a moanin model
#' data(shoemaker2015)
#' data = shoemaker2015$data
#' meta = shoemaker2015$meta
#' moanin_model = create_moanin_model(meta, degrees_of_freedom=6)
#'
#' # The moanin model contains all the information for plotting purposes. The
#' # plot_splines_function will automatically fit the splines from the
#' # information contained in the moanin model
#' plot_splines_data(
#'	data, moanin_model,
#'	subset_data=c("AK050122", "AK043921"),
#'	mfrow=c(2, 2))
#'
#' # The splines can also be smoothed
#' plot_splines_data(data, moanin_model,
#'		     subset_data=c("AK050122", "AK043921"),
#'		     smooth=TRUE, mfrow=c(2, 2))
#' @export
#' @importFrom graphics axis plot.new plot.window
#' @importFrom methods new

plot_splines_data = function(data, moanin_model, colors=NULL, smooth=FALSE,
			     legend=TRUE, legendArgs=NULL, subset_conditions=NULL,
			     subset_data=NULL,
			     simpleY=TRUE,
			     mar=c(2.5, 2.5, 3.0, 1.0),
			     mfrow=NULL, addToPlot=NULL, ...){

    # Start by subsetting the data
    if(!is.null(subset_data)){
	mask = subset_data %in% row.names(data)
	if(!all(mask)){
	    msg = paste0(
		"Some elements of subset_data are not match row.names of data.")
	    stop(msg)
	}
	data = data[subset_data, ]
    }

    n_observations = dim(data)[1]
    n_plots = if(legend) n_observations+1 else n_observations
    if(!is.null(mfrow)){
	if(length(mfrow) != 2){
	    msg = sprintf(
		paste0("Invalid value for argument mfrow. Should be a vector of length 2.",
		       "A vector of length %s was provided.", length(mfrow)))
	    stop(msg)

	}
        if(mfrow[1]*mfrow[2] < n_plots){
	    msg = sprintf(
		paste0(
		    "Invalid value for argument mfrow. Should result in ",
		    "grid for at least %s plots (including a plot for the ",
		    "legend, if legend=TRUE)"),
		n_plots) 
	    stop(msg)
	}
    }else{
        if(n_plots <= 3){ 
            mfrow = c(n_plots, 1)
        }else if(n_plots <= 6){
            mfrow = c(ceiling(n_plots / 2), 2)
        }else if(n_plots <= 12){
            mfrow = c(ceiling(n_plots / 3), 3)
        }else{
            nrow = round(n_plots ** 0.5)
            ncol = ceiling(n_plots / nrow)
            mfrow=c(ncol, nrow)
        }        
    }
    graphics::par(mfrow=mfrow, mar=mar)
    bottomPlots = seq(to=n_observations, by=1, length=mfrow[2])
    sidePlots = seq(from=1, to=n_observations, by=mfrow[2])

   
    ## For legend:
    if(is.null(colors)){
	meta = moanin_model$meta
        groups = levels(meta$Group)
        colors = viridis::viridis(length(groups))
        names(colors) = groups
    }

    plot_names = row.names(data)
    name = NULL 

    # Now plot the different data
    for(i in 1:n_observations){
	if(!is.null(plot_names)){
	    name = plot_names[i]
	}
        plot_centroid_individual(
	    as.vector(data[i, ]),
	    moanin_model, colors=colors,
	    smooth=smooth,
	    subset_conditions=subset_conditions,
	    main=name,
            xaxt=if(!i %in% bottomPlots) "n" else "s",
            yaxt=if(!i %in% sidePlots & simpleY) "n" else "s", ...)
                 
        if(is.function(addToPlot)){
            addToPlot()
        }

        if(!i %in% bottomPlots){
	    axis(1, labels=FALSE)
	}

        if(!i %in% sidePlots & simpleY){
	    axis(2, labels=FALSE)
	}
        
    }

    if(legend){
        plot.new()
        plot.window(
	    xlim=c(0,1),
	    ylim=c(0,1),
	    bty="n",
	    xaxt="n",
	    yaxt="n")

        if(!is.null(subset_conditions)){
	    colors = colors[names(colors) %in% subset_conditions]
	}
        do.call(
	    "legend",
	    c(list(x="center",
		   legend=names(colors),
		   fill=colors, bty="n"),
		   legendArgs))
    }
}

plot_centroid_individual = function(centroid, moanin_model,
			    colors, smooth=FALSE, subset_conditions=NULL, ...){
    meta = moanin_model$meta
    groups = levels(meta$Group)

    if(!is.null(subset_conditions)){
        if(!all(subset_conditions %in% groups)){
	    msg = paste0(
		"subset_conditions argument given by user does not match names",
		" of groups.")
	    stop(msg)
        }else{
	    groups = subset_conditions
	}
    }

    xrange = range(meta$Timepoint)
    if(is.null(dim(centroid))){
        centroid = t(as.matrix(centroid))
    }
    yrange = range(centroid[, meta$Group %in% groups])

    graphics::plot(xrange, yrange, type="n", ...)
    if(smooth){
	# FIXME this is supposed to be on the fitted lines, but I'm not able
	# to get this to work fine in R.
	meta_prediction = create_meta_prediction(moanin_model)
	centroid_fitted = fit_predict_splines(
	    centroid, moanin_model,
	    meta_prediction=meta_prediction)
    }else{
	centroid_fitted = fit_predict_splines(
	    centroid, moanin_model)
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
