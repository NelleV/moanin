setGeneric("plot_splines_data", 
           function(object,data,...) { standardGeneric("plot_splines_data")})

#' Plotting splines
#'
#'@param data matrix containing the data to be plotted, where each row of the
#'  data provided will be plotted as a separate plot. If missing, will rely on
#'  data in \code{assay(object)}
#'@inheritParams DE_timecourse
#'@param colors vector, optional, default NULL. Vector of colors
#'@param smooth boolean, optional, default: FALSE. Whether to smooth the
#'  centroids or not.
#'@param mar vector of margins to set the space around each plot (see
#'  \code{\link{par}})
#'@param legend boolean whether to include a legend (default:TRUE)
#'@param legendArgs list of arguments to be passed to legend command (if
#'  \code{legend=TRUE})
#'@param simpleY boolean, if true, will plot all genes on same y-axis and
#'  minimize the annotation of the y axis to only label the axis in the exterior
#'  plots (the x-axis is always assumed to be the same across all plots and will
#'  always be simplified)
#'@param scale_centroid determines whether the centroid data given in
#'  \code{centroid} should be rescaled to match that of the data
#'  (\code{"toData"}), or the data scaled to match that of centroid
#'  (\code{"toCentroid"}), or simply plotted as is (\code{"none"}).
#'@param subset_conditions list if provided, only plots the subset of conditions
#'  provided. Else, plots all conditions
#'@param centroid numeric vector (or matrix of 1 row) with data to use to fit
#'  the splines. If \code{NULL}, the splines plotted will be from the data.
#'@param subset_data list if provided, only plots the subset of data (ie, the
#'  rows) provided. Can be any valid vector for subsetting a matrix. See
#'  details.
#'@param mfrow a vector of integers of length 2 defining the grid of plots to be
#'  created (see \code{\link{par}}). If missing, the function will set a value.
#'@param addToPlot A function that will be called after the plotting, allowing
#'  the user to add more to the plot.
#'@param xlab label for the x-axis
#'@param ylab label for the y-axis
#' @param xaxis Logical, whether to add x-axis labels to plot (if FALSE can be manually created by user with call to addToPlot)
#' @param yaxis Logical, whether to add y-axis labels to plot (if FALSE can be manually created by user with call to addToPlot)
#'@param ... arguments to be passed to the individual plot commands (Will be
#'  sent to all plot commands)
#' @details If \code{data} is NULL, the data plotted will be from
#'   \code{assay(object)}, after log-transformation if
#'   \code{log_transform(object)=TRUE}.
#' @details If \code{centroid} is missing, then splines will be estimated (per
#'   group) for the the data in \code{data} -- separately for each row of
#'   \code{data}. If \code{centroid} is provided, this data will be used to plot
#'   a spline function, and this same spline will be plotted for each row of
#'   \code{data}. This is useful, for example, in plotting cluster centroids
#'   over a series of genes.
#' @details If the user set \code{log_transform=TRUE} in the creation of the
#'   \code{Moanin} object, the data will be log transformed before plotting and
#'   calculating the spline fits.
#'@return This function creates a plot and does not return anything to the user.
#' @examples
#' # First, load some data and create a moanin model
#' data(exampleData)
#' moanin <- create_moanin_model(data=testData,meta=testMeta, 
#'    degrees_of_freedom=6)
#'
#' # The moanin model contains all the information for plotting purposes. The
#' # plot_splines_data will automatically fit the splines from the
#' # information contained in the moanin model
#' genes <- c("NM_001042489", "NM_008725")
#' plot_splines_data(moanin, subset_data=genes,
#' mfrow=c(2, 2))
#' # By default, same axis for all genes. Can change with 'simpleY=FALSE'
#' plot_splines_data(moanin, subset_data=genes,
#'    smooth=TRUE, mfrow=c(2,2), simpleY=FALSE)   
#'
#' # The splines can also be smoothed
#' plot_splines_data(moanin, subset_data=genes,
#'    smooth=TRUE, mfrow=c(2, 2))
#' # You can provide different data (on same subjects),
#' # instead of data in moanin object
#' # (in which case moanin just provides grouping information)
#' plot_splines_data(moanin, data=1/assay(moanin), subset_data=genes,
#'    smooth=TRUE, mfrow=c(2, 2))
#'    
#' # You can also provide data to use for fitting splines to argument  
#' # "centroid". This is helpful for overlaying centroids or predicted data
#' # Here we do a silly example, just to demonstrate syntax, 
#' # where we use the data from the first gene as our centroid to fit a
#' # spline estimate, but plot data from genes 3-4
#' plot_splines_data(moanin, centroid=assay(moanin[1,]), subset_data=3:4,
#'    smooth=TRUE, mfrow=c(2,2))
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
                simpleY=TRUE, centroid=NULL,
                scale_centroid=c("toData","toCentroid","none"),
                mar=c(2.5, 2.5, 3.0, 1.0),
                mfrow=NULL, addToPlot=NULL, ylab="", xaxis=TRUE,
                yaxis=TRUE, xlab="Time",...){
    scale_centroid <- match.arg(scale_centroid)
    check_data_meta(data=data,object=object)
    if(!is.null(centroid)){
        check_data_meta(data=centroid,object=object)
        if(!is.null(dim(centroid)) && nrow(centroid)>1){
            stop("Centroid must be a single vector (or matrix of 1 row) to data, for fitting the functional form")
        } 
    }
    if(!is.null(subset_data) )  data <- data[subset_data,]
    
    n_observations <- dim(data)[1]
    ### Work out the mfrow/number of plots and check makes sense
    n_plots <- if(legend) n_observations+1 else n_observations
    if(!is.null(mfrow)){
        if(length(mfrow) != 2){
            msg <- sprintf(
                paste0("Invalid value for argument mfrow.",
                       "Should be a vector of length 2.",
                       "A vector of length %s was provided.", length(mfrow)))
            stop(msg)
            
        }
        if(mfrow[1]*mfrow[2] < n_plots){
            msg <- sprintf(
                paste0(
                    "Invalid value for argument mfrow. Should result in ",
                    "grid for at least %s plots (including a plot for the ",
                    "legend, if legend=TRUE)"),
                n_plots) 
            stop(msg)
        }
    }else{
        if(n_plots <= 3){ 
            mfrow <- c(n_plots, 1)
        }else if(n_plots <= 6){
            mfrow <- c(ceiling(n_plots / 2), 2)
        }else if(n_plots <= 12){
            mfrow <- c(ceiling(n_plots / 3), 3)
        }else{
            nrow <- round(n_plots ** 0.5)
            ncol <- ceiling(n_plots / nrow)
            mfrow=c(ncol, nrow)
        }        
    }
    graphics::par(mfrow=mfrow, mar=mar)
    bottomPlots <- seq(to=n_observations, by=1, length=mfrow[2])
    sidePlots <- seq(from=1, to=n_observations, by=mfrow[2])
    
    
    ## For legend:
    if(is.null(colors)){
        groups <- levels(group_variable(object))
        colors <- viridis::viridis(length(groups))
        names(colors) <- groups
    }
    
    plot_names <- row.names(data)
    name <- NULL 
    if(!is.null(centroid)){
        if(scale_centroid=="toData"){
            centroid <- align_data_onto_centroid(data, centroid, 
                                                returnType="centroid")
        }
        if(scale_centroid=="toCentroid"){
            data <- align_data_onto_centroid(data, centroid, returnType="data")
        }
        if(scale_centroid %in% c("toCentroid","none")){
            #make it a matrix
            centroid <- matrix(centroid, nrow=nrow(data),ncol=length(centroid),
                              byrow=TRUE)
        }
    }
    
    ## Fix up the y-axis range
    if("ylim" %in% names(list(...))){
        ylim <- list(...)$ylim
    }
    else ylim <- NULL
    if(simpleY & is.null(ylim)){
        # use same y-axis range
        yrange <- range(data)
    }
    else{
        if(!is.null(ylim)) yrange <- ylim
        else yrange <- NULL
    }

    # Now plot the different data
    for(i in seq_len(n_observations)){
        if(!is.null(plot_names)){
            name <- plot_names[i]
        }
        xaxt<-if(!i %in% bottomPlots & xaxis) "n" else "s"
        if("xaxt" %in% names(list(...))) xaxt<-list(...)$xaxt
        yaxt<-if(!i %in% bottomPlots & yaxis) "n" else "s"
        if("yaxt" %in% names(list(...))) yaxt<-list(...)$yaxt
        plot_centroid_individual(
            data=data[i, ], centroid=if(is.null(centroid)) centroid else centroid[i,],
            object[i,], colors=colors,
            smooth=smooth,
            subset_conditions=subset_conditions,
            main=name, yrange=yrange,
            xaxt=xaxt,
            yaxt=yaxt, 
            xlab=xlab,ylab=ylab,...)
        
        if(is.function(addToPlot)){
            addToPlot()
        }
        
        if(!i %in% bottomPlots & xaxis){
            axis(1, labels=FALSE)
        }
        
        if(!i %in% sidePlots & simpleY & yaxis){
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
            colors <- colors[names(colors) %in% subset_conditions]
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

#' @aliases plot_splines_data,Moanin,DataFrame-method
#' @rdname plot_splines_data
#' @export
setMethod("plot_splines_data",c("Moanin","DataFrame"),
          function(object, data, ...){
              plot_splines_data(object,data=data.matrix(data),...)
          }
)
#' @aliases plot_splines_data,Moanin,missing-method
#' @rdname plot_splines_data
#' @export
setMethod("plot_splines_data",c("Moanin","missing"),
        function(object, data, ...){
            data<-get_log_data(object)
            plot_splines_data(object,data=data,...)
        }
)

# centroid and data are vectors
# moanin_model will separate them into groups...
plot_centroid_individual <- function(data, centroid, moanin_model,
                                    colors, smooth=FALSE, 
                                    subset_conditions=NULL,yrange=NULL,...){
    if(is.null(data)){
        data=as.vector(get_log_data(moanin_model))
    }
    if(is.null(centroid)) centroid<-data
    if(!inherits(moanin_model,"Moanin")) 
        stop("Internal coding error: expecting Moanin class")
    gpVar <- group_variable_name(moanin_model)
    tpVar <- time_variable_name(moanin_model)
    groups <- levels(group_variable(moanin_model))
    
    if(!is.null(subset_conditions)){
        if(!all(subset_conditions %in% groups)){
            msg <- paste0(
                "subset_conditions argument given by user does not match names",
                " of groups.")
            stop(msg)
        }else{
            groups <- subset_conditions
        }
    }
    
    xrange <- range(time_variable(moanin_model))
    if(inherits(centroid,"DataFrame")) centroid<-as.matrix(centroid)
    if(is.null(dim(centroid))) centroid <- t(as.matrix(centroid))
    if(is.null(yrange)){
        yrangeCentroid <- range(centroid[, group_variable(moanin_model) %in% groups])
        yrangeData <- range(data[group_variable(moanin_model) %in% groups])
        yrange <- range(c(yrangeCentroid,yrangeData))
    }
    graphics::plot(xrange, yrange, type="n", ...)
    if(smooth){
        # FIXME this is supposed to be on the fitted lines, but I'm not able
        # to get this to work fine in R.
        meta_prediction <- create_meta_prediction(moanin_model)
        centroid_fitted <- fit_predict_splines(
            centroid, moanin_model,
            meta_prediction=meta_prediction)
    }else{
        centroid_fitted <- fit_predict_splines(
            centroid, moanin_model)
    }
    
    
    # scatter points for values
    for(i in seq_along(groups)){
        group <- groups[i]
        color <- colors[group]
        
        # Start by individual points
        mask <- group_variable(moanin_model) == group
        time <- time_variable(moanin_model)[mask]
        indx <- order(time)
        graphics::lines(time[indx], data[mask][indx], type="p",
                        col=color, pch=16,
                        lwd=0)
        if(smooth){
            mask <- meta_prediction[,gpVar] == group
            time <- meta_prediction[,tpVar][mask]
            indx <- order(time)
        }
        graphics::lines(time[indx], centroid_fitted[mask][indx], type="l",
                        col=color, lwd=1)
        
    }
}
