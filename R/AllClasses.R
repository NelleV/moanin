setOldClass("formula", S4Class="formula")
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
#' important methods have their own help page. Simple helper methods are described in the
#' Methods section below. For a comprehensive list of methods specific to this class
#' see the Reference Manual.
#' @slot time_variable_name character value giving the column in \code{colData} that defines the time variable (must be of class \code{numeric})
#' @slot group_variable_name character value giving the column in \code{colData} that defines the grouping variable (must be of class \code{factor})
#' @slot basis A basis matrix, where each row
#'  corresponds to the evaluation of a sample on the basis function (thus one
#'  column for each basis function).
#' @slot formula a formula. The formula used in creating the basis matrix
#' @slot degrees_of_freedom a numeric integer. Number of degrees of freedom used in creating basis matrix. If NULL, degrees of freedom is not known (usually if user provided basis without degrees of freedom)
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
        formula="formulaOrNULL", 
        basis="matrix",
        group_variable_name="character",
        time_variable_name="character",
        degrees_of_freedom="numericOrNULL"
    )
)

setValidity("Moanin", function(object) {
    if(! group_variable_name(object) %in% colnames(colData(object))) return("group_variable_name slot must match the name of one of the columns of colData")
    else{
        if(!is.factor(group_variable(object))) return(group_variable_name(object),"is not a numeric column in the colData of the object")
    }
    if(! time_variable_name(object) %in% colnames(colData(object))) return("time_variable_name slot must match the name of one of the columns of colData")
    else{
        if(!is.numeric(time_variable(object))) return(time_variable_name(object),"is not a numeric column in the colData of the object")
    }

    return(TRUE)
})

#' Create a Moanin object
#'@description The constructor \code{create_moanin_model} creates an object of the
#'  class \code{Moanin}. 
#'
#' 
#' @param group_variable_name A character value giving the column that corresponds to
#'   the grouping variable to test for DE. By default "Group"
#' @param time_variable_name A character value giving the column that corresponds to
#'   the time variable. By default "Timepoint".
#' @param formula formula object, optional, default: NUlL. Used to construct
#'   splines from the data in \code{meta}. See details.
#'@param basis	matrix, optional, default: NULL. A basis matrix, where each row
#'  corresponds to the evaluation of a sample on the basis function (thus one
#'  column for each basis function).
#'@param degrees_of_freedom int, optional. Number of degrees of freedom to use
#'  if neither the basis nor the formula is provided. If not provided by the
#'  user, internally will be set to 4
#'@details If neither \code{formula} nor \code{basis} is given, then by default,
#'  the function will create a basis matrix based on the formula:
#'  \preformatted{formula = ~Group:ns(Timepoint, df=4) + Group + 0}
#'@details Note that the meta data will have levels dropped (via \code{droplevels}). 
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
#' moanin = create_moanin_model(data=testData,meta=testMeta, degrees_of_freedom=6)
#' moanin
#' @export
#' @importFrom splines ns
#' @importFrom stats as.formula
#' @rdname Moanin-class
setGeneric(
    name = "create_moanin_model",
    def = function(data,  ...) {
        standardGeneric("create_moanin_model")
    }
)
#' @rdname Moanin-class
#' @param meta if \code{data} is of class matrix, argument \code{meta} must be given a \code{data.frame} containing the metadata in columns, and rows. This will be made into the \code{colData} of the resulting \code{Moanin} object.
#'   corresponding to different samples. 
#' @export
setMethod(
    f = "create_moanin_model",
    signature = signature("matrix"),
    definition = function(data, meta, ...){
        create_moanin_model(SummarizedExperiment(data, colData=meta,...))
    })
#' @rdname ClusterExperiment-class
setMethod(
    f = "create_moanin_model",
    signature = signature("SummarizedExperiment"),
    definition = function(data, formula=NULL, basis=NULL,
                               group_variable_name="Group",time_variable_name="Timepoint",
                               degrees_of_freedom=NULL){

    if(!is.null(basis) & !is.null(formula)){
        msg = paste("both basis and formula ",
                    "are provided by the user. Please provide one or ",
                    "the other", sep="")
        stop(msg)
    }
    
    if(is.null(basis)){
        if(is.null(formula)){
            if(is.null(degrees_of_freedom)){
                degrees_of_freedom = 4
            }
            formulaText<-paste0("~",group_variable_name," + ",group_variable_name,":splines::ns(",time_variable_name,",df=",degrees_of_freedom,") + 0")
            formula = stats::as.formula(formulaText)# (
            #                 ~Group + Group:splines::ns(Timepoint, df=degrees_of_freedom) + 0)
        }
        basis = stats::model.matrix(formula, data=meta)
    }else{
        basis = as.matrix(basis)
    }
    
    splines_model = new("Moanin", 
                        degrees_of_freedom=degrees_of_freedom,
                        basis=basis,
                        formula=formula,
                        time_variable_name=time_variable_name,
                        group_variable_name=group_variable_name
                        )
    group_variable(splines_model)<-droplevels(group_variable(splines_model))
    # Just create this one.
    if(!("WeeklyGroup" %in% colnames(colData(splines_model)))){
        colData(splines_model)$WeeklyGroup = as.factor(
            make.names(group_variable(splines_model):as.factor(time_variable(splines_model))))
    }
    #somewhat wasteful here, because doing checks twice!
    validObject(splines_model)
    return(splines_model)
    }
)


