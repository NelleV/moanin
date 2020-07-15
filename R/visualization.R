setGeneric("plot_splines_data", 
           function(object,data,...) { standardGeneric("plot_splines_data")})

#' Plotting splines
#'
#'@param data matrix containing the data to be plotted (such as centroids or
#'  data), where each row of the data provided will be plotted as a separate
#'  plot. If missing, will rely on data in \code{assay(object)}
#'@inheritParams DE_timecourse
#'@param colors vector, optional, default NULL. Vector of colors
#'@param smooth boolean, optional, default: FALSE. Whether to smooth the
#'  centroids or not.
#'@param mar vector of margins to set the space around each plot (see
#'  \code{\link{par}})
#'@param legend boolean whether to include a legend (default:TRUE)
#'@param legendArgs list of arguments to be passed to legend command (if
#'  \code{legend=TRUE})
#'@param simpleY boolean, if true, will minimize the annotation of the y axis to
#'  only label the axis in the exterior plots (the x-axis is always assumed to
#'  be the same across all plots and will always be simplified)
#'@param subset_conditions list if provided, only plots the subset of conditions
#'  provided. Else, plots all conditions
#'@param data_smooth_only logical. If TRUE, then the data given in \code{data}
#'  will only be used to plot the functional form, and the data from the
#'  \code{Moanin} object will be used for the points. See details.
#'@param subset_data list if provided, only plots the subset of data (ie, the
#'  rows) provided. Can be any valid vector for subsetting a matrix. See details.
#'@param mfrow a vector of integers of length 2 defining the grid of plots to be
#'  created (see \code{\link{par}}). If missing, the function will set a value.
#'@param addToPlot A function that will be called after the plotting, allowing
#'  the user to add more to the plot.
#'@param ... arguments to be passed to the individual plot commands (Will be
#'  sent to all plot commands)
#' @details If \code{data_smooth_only=TRUE}, then the data plotted will be from
#'   \code{assay(object)} and the matrix of \code{data} will only be used to fit
#'   functional forms. This is useful, for example, in plotting cluster
#'   centroids over the data in \code{object}. Note that in this case,
#'   \code{subset_data} will be applied to \code{object} and NOT to \code{data}.
#'   This allows you to plot the same centroids over a series of genes.
#'@return This function creates a plot and does not return anything to the user.
#' @examples
#' # First, load some data and create a moanin model
#' data(exampleData)
#' moanin = create_moanin_model(data=testData,meta=testMeta, 
#'    degrees_of_freedom=6)
#'
#' # The moanin model contains all the information for plotting purposes. The
#' # plot_splines_function will automatically fit the splines from the
#' # information contained in the moanin model
#' genes = c("NM_001042489", "NM_008725")
#' plot_splines_data(moanin, subset_data=genes,
#' mfrow=c(2, 2))
#'
#' # The splines can also be smoothed
#' plot_splines_data(moanin, subset_data=genes,
#'    smooth=TRUE, mfrow=c(2, 2))
#'    
#' # You can provide different data on same subjects,
#' # instead of data in moanin object
#' # (in which case moanin just provides grouping information)
#' plot_splines_data(moanin, data=1/(assay(moanin)), subset_data=genes,
#'    smooth=TRUE, mfrow=c(2, 2))
#'    
#' # You can also provide different data for fitting splines, 
#' # but plot data in Moanin object
#' # This is helpful for overlaying centroids or predicted data
#' # Here we do a silly example, where we use the data from genes 1-2 to fit 
#' # splines, but plot data from genes 3-4, just to demonstrate syntax
#' plot_splines_data(moanin, data=assay(moanin[1:2,]), subset_data=3:4,
#'    smooth=TRUE, mfrow=c(2, 2), data_smooth_only=TRUE)
#' @export
#' @name plot_splines_data
#' @aliases plot_splines_data,Moanin,matrix-method
#' @importFrom graphics axis plot.new plot.window
#' @importFrom methods new
#' @importFrom viridis viridis
setMethod("plot_splines_data",c("Moanin","matrix"),
    function(object, data, colors=NULL, smooth=FALSE,
                legend=TRUE, legendArgs=NULL, 
                subset_conditions=NULL,
                subset_data=NULL,
                simpleY=TRUE, data_smooth_only=FALSE,
                mar=c(2.5, 2.5, 3.0, 1.0),
                mfrow=NULL, addToPlot=NULL, ...){
    
    check_data_meta(data=data,object=object)
    if(!is.null(subset_data) ){
        if(data_smooth_only){
            object = object[subset_data,]
        }
        else{
            data = data[subset_data,]
        }
    }
    n_observations = if(data_smooth_only) dim(object)[1] else dim(data)[1]
    ### Work out the mfrow/number of plots and check makes sense
    n_plots = if(legend) n_observations+1 else n_observations
    if(!is.null(mfrow)){
        if(length(mfrow) != 2){
            msg = sprintf(
                paste0("Invalid value for argument mfrow.",
                       "Should be a vector of length 2.",
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
        groups = levels(group_variable(object))
        colors = viridis::viridis(length(groups))
        names(colors) = groups
    }
    
    plot_names = row.names(data)
    name = NULL 
    
    # Now plot the different data
    for(i in seq_len(n_observations)){
        if(!is.null(plot_names)){
            name = plot_names[i]
        }
        plot_centroid_individual(
            centroid=as.vector(data[i, ]),
            object[i,], colors=colors,
            smooth=smooth,
            subset_conditions=subset_conditions,
            main=name, data=if(data_smooth_only) NULL else as.vector(data[i, ]),
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
)

#' @aliases plot_splines_data,Moanin,numeric-method
#' @rdname plot_splines_data
#' @export
setMethod("plot_splines_data",c("Moanin","numeric"),
          function(object, data, ...){
              plot_splines_data(object,data=matrix(data,nrow=1),...)
          }
)
#' @aliases plot_splines_data,Moanin,data.frame-method
#' @rdname plot_splines_data
#' @export
setMethod("plot_splines_data",c("Moanin","data.frame"),
          function(object, data, ...){
              plot_splines_data(object,data=data.matrix(data),...)
          }
)
#' @aliases plot_splines_data,Moanin,missing-method
#' @rdname plot_splines_data
#' @export
setMethod("plot_splines_data",c("Moanin","missing"),
        function(object, data, ...){
            data<-assay(object)
            if(object@log_transform) log(data+1)
              plot_splines_data(object,data=data,data_smooth_only=FALSE,...)
        }
)

# centroid and data are vectors
# moanin_model will separate them into groups...
plot_centroid_individual = function(centroid, moanin_model,
                                    colors, smooth=FALSE, 
                                    subset_conditions=NULL,
                                    data_smooth_only, data,...){
    if(is.null(data)){
        if(moanin_model@log_transform) data = as.vector(log(assay(moanin_model)+1))
        else data = as.vector(assay(moanin_model) )
    }
    if(!inherits(moanin_model,"Moanin")) 
        stop("Internal coding error: expecting Moanin class")
    gpVar = group_variable_name(moanin_model)
    tpVar = time_variable_name(moanin_model)
    groups = levels(group_variable(moanin_model))
    
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
    
    xrange = range(time_variable(moanin_model))
    if(is.null(dim(centroid))){
        centroid = t(as.matrix(centroid))
    }
    yrangeCentroid = range(centroid[, group_variable(moanin_model) %in% groups])
    yrangeData = range(data[group_variable(moanin_model) %in% groups])
    yrange = range(c(yrangeCentroid,yrangeData))
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
    for(i in seq_along(groups)){
        group = groups[i]
        color = colors[group]
        
        # Start by individual points
        mask = group_variable(moanin_model) == group
        time = time_variable(moanin_model)[mask]
        indx = order(time)
        graphics::lines(time[indx], data[mask][indx], type="p",
                        col=color, pch=16,
                        lwd=0)
        if(smooth){
            mask = meta_prediction[,gpVar] == group
            time = meta_prediction[,tpVar][mask]
            indx = order(time)
        }
        graphics::lines(time[indx], centroid_fitted[mask][indx], type="l",
                        col=color, lwd=1)
        
    }
} 

# DELETE ME
# plot_gene_splines = function(data, meta, gene_name, colors=NULL){
#     # First, select the gene:
#     if(!(gene_name %in% row.names(data))){
#         msg = paste("moanin::plot_gene_splines: The gene_name provided '",
#                     gene_name, "' is not in the data", sep="")
#     }
#     gene_data = data[gene_name, ]
# }
