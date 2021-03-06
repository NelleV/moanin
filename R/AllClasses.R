#setOldClass("formula", S4Class="formula")
setClassUnion("numericOrNULL",members=c("numeric", "NULL"))
setClassUnion("formulaOrNULL",members=c("formula", "NULL"))

#' @title Class Moanin
#'
#' @description \code{Moanin} is a class that extends
#'   \code{SummarizedExperiment} and is used to store the additional spline
#'   basis and meta data for timecourse analysis.
#'
#' @docType class
#' @aliases Moanin Moanin-class
#' @description In addition to the slots of the \code{SummarizedExperiment}
#' class, the \code{Moanin} object has the additional slots described
#' in the Slots section.
#'
#' @description There are several methods implemented for this class. The most
#'   important methods have their own help page. Simple helper methods are
#'   described in the Methods section below. For a comprehensive list of methods
#'   specific to this class see the Reference Manual.
#' @slot time_variable_name character value giving the column in \code{colData}
#'   that defines the time variable (must be of class \code{numeric})
#' @slot group_variable_name character value giving the column in \code{colData}
#'   that defines the grouping variable (must be of class \code{factor})
#' @slot basis_matrix A basis matrix, where each row corresponds to the
#'   evaluation of a sample on the basis function (thus one column for each
#'   basis function).
#' @slot spline_formula a formula. The formula used in creating the basis matrix
#' @slot degrees_of_freedom a numeric integer. Number of degrees of freedom used
#'   in creating basis matrix. If NULL, degrees of freedom is not known (usually
#'   if user provided basis without degrees of freedom)
#' @slot log_transform logical, whether to log-transform the data for certain
#'   methods
#' @name Moanin-class
#' @aliases Moanin
#' @rdname Moanin-class
#' @import SummarizedExperiment
#' @import methods
#' @export
#'
setClass(
    Class = "Moanin",
    contains = "SummarizedExperiment",
    slots = list(
        spline_formula="formulaOrNULL", 
        basis_matrix="matrix",
        group_variable_name="character",
        time_variable_name="character",
        degrees_of_freedom="numericOrNULL",
        log_transform="logical"
    )
)

setValidity("Moanin", function(object) {
    if(! group_variable_name(object) %in% colnames(colData(object))) 
        return(paste("group_variable_name slot must match the name",
            "of one of the columns of colData"))
    else{
        if(!is.factor(group_variable(object))) 
            return(paste(group_variable_name(object),
                "is not a numeric column in the colData of the object"))
    }
    if(! time_variable_name(object) %in% colnames(colData(object))) 
        return(paste("time_variable_name slot must match the", 
            "name of one of the columns of colData"))
    else{
        if(!is.numeric(time_variable(object))) 
            return(paste(time_variable_name(object),
                "is not a numeric column in the colData of the object"))
    }
    if(nrow(basis_matrix(object)) != ncol(object)) 
        return(paste("number of rows of the basis_matrix doesn't match",
            "number of samples of object"))

    return(TRUE)
})

setGeneric(
    name = "create_moanin_model",
    def = function(data,  ...) {
        standardGeneric("create_moanin_model")
    }
)
#'Create a Moanin object
#'@description The constructor \code{create_moanin_model} creates an object of
#'  the class \code{Moanin}.
#'
#'@param data The input data. Can be a \code{\link{SummarizedExperiment}} class,
#'  or matrix/data.frame. If the input data is a \code{matrix} or
#'  \code{data.frame}, then the user must also provide input to the \code{meta}
#'  argument, which will be transformed into \code{colData} of the resulting
#'  \code{Moanin} object
#'@param meta Meta data on the samples (columns) of the \code{data} argument.
#'  Must be given f input \code{data} is a \code{matrix} or \code{data.frame}.
#'  If input is \code{SummarizedExperiment}, this argument is ignored.
#'@param group_variable_name A character value giving the column that
#'  corresponds to the grouping variable to test for DE. By default "Group"
#'@param time_variable_name A character value giving the column that corresponds
#'  to the time variable. By default "Timepoint".
#'@param spline_formula formula object, optional, default: NUlL. Used to
#'  construct splines from the data in \code{meta}. See details.
#'@param basis_matrix matrix, optional, default: NULL. A basis matrix, where
#'  each row corresponds to the evaluation of a sample on the basis function
#'  (thus one column for each basis function).
#'@param degrees_of_freedom int, optional. Number of degrees of freedom to use
#'  if neither the basis_matrix nor the spline_formula is provided. If not
#'  provided by the user, internally will be set to 4
#'@param log_transform whether the data should be log-transformed by certain
#'  methods (see \code{\link{splines_kmeans}})
#'@param drop_levels Logical, whether to perform \code{\link{droplevels}} on the
#'  grouping variable (i.e. remove empty levels)
#'@param ... arguments passed from methods to the \code{SummarizedExperiment}
#'  method.
#'@details If neither \code{spline_formula} nor \code{basis_matrix} is given,
#'  then by default, the function will create a basis matrix based on the
#'  formula: \preformatted{spline_formula = ~Group:ns(Timepoint, df=4) + Group +
#'  0}
#'@details Note that the meta data will have levels dropped (via 
#'\code{droplevels}). 
#'@details Input to \code{data} that is given as a class \code{DataFrame} or
#'  \code{data.frame} will be converted to class \code{matrix}. The reason
#'  for this is that use of a \code{data.frame} creates errors in taking duplicate rows/columns of \code{SummarizedExperiment}, as in bootstrapping. 
#'  Users who absolutely want the object to hold a
#'  object that is not a matrix can construct a \code{SummarizedExperiment} object
#'  (which will not convert the input into a \code{matrix}), and
#'  use this as input to \code{create_moanin_model}.
#'@return An object of class \code{Moanin}
#' @examples
#' # Load some data
#' data(exampleData)
#'
#' # Use the default options
#' moanin = create_moanin_model(data=testData,meta=testMeta)
#' moanin
#'
#' # Change the number of degrees of freedom
#' moanin = create_moanin_model(data=testData,meta=testMeta, 
#'    degrees_of_freedom=6)
#' moanin
#' @export
#' @importFrom splines ns
#' @importFrom stats as.formula
#' @rdname Moanin-class
#' @aliases create_moanin_model
#' @aliases Moanin
#' @importClassesFrom S4Vectors DataFrame
#' @export
setMethod(
    f = "create_moanin_model",
    signature = signature("DataFrame"),
    definition = function(data, meta, ...){
        create_moanin_model(data.matrix(data),...)
    })
#' @rdname Moanin-class
#' @export
setMethod(
    f = "create_moanin_model",
    signature = signature("data.frame"),
    definition = function(data, ...){
        create_moanin_model(data.matrix(data), ...)
    })
#' @rdname Moanin-class
#' @export
setMethod(
    f = "create_moanin_model",
    signature = signature("matrix"),
    definition = function(data,  meta,...){
        if(missing(meta)) stop("Must provide argument meta if input to data is",
                "matrix-like")
        create_moanin_model(SummarizedExperiment(data, colData=meta),...)
    })
#' @rdname Moanin-class
#' @export
setMethod(
    f = "create_moanin_model",
    signature = signature("SummarizedExperiment"),
    definition = function(data, spline_formula=NULL, basis_matrix=NULL,
                group_variable_name="Group",time_variable_name="Timepoint",
                degrees_of_freedom=NULL,log_transform=FALSE,drop_levels=TRUE){

    if(!is.null(basis_matrix) & !is.null(spline_formula)){
        msg = paste("both basis_matrix and spline_formula ",
                    "are provided by the user. Please provide one or ",
                    "the other", sep="")
        stop(msg)
    }
    if(!("WeeklyGroup" %in% colnames(colData(data)))){
        colData(data)$WeeklyGroup = as.factor(
            make.names(colData(data)[,group_variable_name]:as.factor(
                colData(data)[,time_variable_name])))
    }
    # Must be done *before* build basis
    if(drop_levels){
        colData(data)[,group_variable_name]<-
            droplevels(colData(data)[,group_variable_name])
    }        
    if(is.null(basis_matrix)){
        if(is.null(spline_formula)){
            if(is.null(degrees_of_freedom)){
                degrees_of_freedom = 4
            }
            formulaText<-paste0("~",group_variable_name," + ",
                                group_variable_name,
                                ":splines::ns(",time_variable_name,
                                ",df=",degrees_of_freedom,") + 0")
            spline_formula = stats::as.formula(formulaText)
        }
        basis_matrix = stats::model.matrix(spline_formula, data=colData(data))
    }else{
        basis_matrix = as.matrix(basis_matrix)
    }
    
    splines_model = new("Moanin", data,
                        degrees_of_freedom=degrees_of_freedom,
                        basis_matrix=basis_matrix,
                        spline_formula=spline_formula,
                        time_variable_name=time_variable_name,
                        group_variable_name=group_variable_name,
                        log_transform=log_transform
                        )
    # Just create this one.
    
    return(splines_model)
    }
)


