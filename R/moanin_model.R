# # Creates a "class" containing the information related to the splines used.
#
# #' Depricated function for creating a moanin_model object
# #'
# #' This function is depricated. Users should use \code{\link{create_moanin_model}}
# #'
# #'@param inheritParams create_moanin_model
# #'@details This function is depricated. Users should use \code{\link{create_moanin_model}}
# #' @export
# #' @keywords internal
# create_splines_model = function(meta, formula=NULL, basis=NULL,
#                 degrees_of_freedom=4){
#     .Deprecated(new="create_moanin_model")
#     return(create_moanin_model(meta, formula=formula, basis=basis,
#                    degrees_of_freedom=degrees_of_freedom))
# }

#' Create a moanin model
#' 
#' This model is used to conserve information relating to the model, such as
#' the formula, the basis, the meta data.
#'
#' @param meta	\code{data.frame} containing the metadata in columns, and rows
#'   corresponding to different samples. 
#' @param group_variable A character value giving the column that corresponds to
#'   the grouping variable to test for DE. By default "Group"
#' @param time_variable A character value giving the column that corresponds to
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
#'@return An object of class \code{moanin_model}; a list with the following elements:
#'\itemize{
#'\item{\code{basis}}{A basis matrix, where each row
#'  corresponds to the evaluation of a sample on the basis function (thus one
#'  column for each basis function). }
#'\item{\code{meta}}{The data frame provided by user to \code{meta}}
#'\item{\code{degrees_of_freedom}}{The degrees of freedom provided by the user 
#' (or set by the function internally)}
#'\item{\code{formula}}{The formula used to construct the basis matrix 
#'(unless basis matrix is provided by the user).}
#'} 
#' @examples
#' # Load some data
#' data(exampleData)
#'
#' # Use the default options
#' moanin = create_moanin_model(testMeta)
#' print(moanin)
#'
#' # Change the number of degrees of freedom
#' moanin = create_moanin_model(testMeta, degrees_of_freedom=6)
#' print(moanin)
#' @export
#' @aliases print.moanin_model
#' @importFrom splines ns
#' @importFrom stats as.formula
create_moanin_model = function(meta, formula=NULL, basis=NULL,
                               group_variable="Group",time_variable="Timepoint",
                               degrees_of_freedom=NULL){
    meta = check_meta(meta,
                      group_variable=group_variable,
                      time_variable=time_variable)
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
            formulaText<-paste0("~",group_variable," + ",group_variable,":splines::ns(",time_variable,",df=",degrees_of_freedom,") + 0")
            formula = stats::as.formula(formulaText)# (
            #                 ~Group + Group:splines::ns(Timepoint, df=degrees_of_freedom) + 0)
        }
        basis = stats::model.matrix(formula, data=meta)
    }else{
        basis = as.matrix(basis)
    }
    
    splines_model = list()
    splines_model$degrees_of_freedom = degrees_of_freedom
    splines_model$"basis" = basis
    splines_model$"meta" = meta
    splines_model$"formula" = formula
    splines_model$time_variable = time_variable
    splines_model$group_variable = group_variable
    class(splines_model) = "moanin_model"
    
    return(splines_model)
}

#' @rdname create_moanin_model
#' @param x a \code{moanin_model} object
#' @param ... arguments passed to generic \code{print}
#' @export
print.moanin_model<-function(x,...){
    N<-nrow(x$meta)
    cat("moanin_model object on",N,"samples containing the following information:\n")
    cat("1) Meta data with",ncol(x$meta),"variables\n")
    if(ncol(x$meta)<=10) print(colnames(x$meta))
    cat(paste0("2) Group variable given by '",x$group_variable,"' with the following levels:\n"))
    print(summary(x$meta$Group))
    cat(paste0("3) Time variable given by '",x$time_variable,"'\n"))
    cat("4) Basis matrix with",ncol(x$basis),"basis functions\n")
    if(!is.null(x$formula)){
        cat("Basis matrix was constructed with the following formula\n")
        form<-gsub("\\s{2,}","",deparse(x$formula)) #get rid of extra spaces
        cat(paste(form,collapse="",sep=""),"\n")
    }
    else{
        if(is.null(x$degrees_of_freedom))
            cat("Basis matrix was provided by user, formula and degrees_of_freedom=NULL\n")
        else
            cat("Basis matrix and degrees of freedom provided by user, equal to",x$degrees_of_freedom,"\n")
    }
    cat("To access this information use:\n")
    cat(paste("\t<object_name>$",names(x),collapse="\n",sep=""))
}
