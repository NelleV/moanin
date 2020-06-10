# Creates a "class" containing the information related to the splines used.

#' Depricated function for creating a moanin_model object
#'
#' This function is depricated. Users should use \code{\link{create_moanin_model}}
#'
#'@param inheritParams create_moanin_model
#'@details This function is depricated. Users should use \code{\link{create_moanin_model}}
#' @export
#' @keywords internal
create_splines_model = function(meta, formula=NULL, basis=NULL,
				degrees_of_freedom=4){
    .Deprecated(new="create_moanin_model")
    return(create_moanin_model(meta, formula=formula, basis=basis,
			       degrees_of_freedom=degrees_of_freedom))
}

#' Create a moanin model
#'
#' This model is used to conserve information relating to the model, such as
#' the formula, the basis, the meta data.
#'
#' @param meta	\code{data.frame} containing the metadata in columns, and rows
#'   corresponding to different samples. Metadata needs to contain the columns
#'   "Group" and "Timepoint".
#' @param formula formula object, optional, default: NUlL. Used to construct splines from the data in \code{meta}. See details.
#'@param basis	matrix, optional, default: NULL. A basis matrix, where each row
#'  corresponds to the evaluation of a sample on the basis function (thus one
#'  column for each basis function).
#'@param degrees_of_freedom int, optional, default: 4. Number of degrees of
#'  freedom to use if neither the basis nor the formula is provided
#'@details If neither \code{formula} nor \code{basis} is given, then by default,
#'  the function will create a basis matrix based on the formula:
#'  \preformatted{formula = ~Group:ns(Timepoint, df=4) + Group + 0}
#'@return An object of class \code{moanin_model}; a list with the following elements:
#'\itemize{
#'\item{\code{basis}}{A basis matrix, where each row
#'  corresponds to the evaluation of a sample on the basis function (thus one
#'  column for each basis function). }
#'\item{\code{meta}}{The data frame provided by user to \code{meta}}
#'\item{\code{degrees_of_freedom}}{The degrees of freedom provided by the user}
#'\item{\code{formula}}{The formula used to construct the basis matrix 
#'(unless basis matrix is provided by the user).}
#'} 
#' @examples
#' # Load some data
#' data(shoemaker2015)
#' meta = shoemaker2015$meta
#'
#' # Use the default options
#' moanin = create_moanin_model(meta)
#' print(dim(moanin$basis))
#'
#' # Change the number of degrees of freedom
#' moanin = create_moanin_model(meta, degrees_of_freedom=6)
#' print(dim(moanin$basis))
#' @export
create_moanin_model = function(meta, formula=NULL, basis=NULL,
                               degrees_of_freedom=4){
    meta = check_meta(meta)
    if(!is.null(basis) & !is.null(formula)){
        msg = paste("moanin::create_splines_model: both basis and formula ",
                    "are provided by the user. Please provide one or ",
                    "the other", sep="")
        stop(msg)
    }
    
    if(is.null(basis)){
        if(is.null(formula)){
            formula = (
                ~Group + Group:splines::ns(Timepoint, df=degrees_of_freedom) + 0)
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
    class(splines_model) = "moanin_model"
    
    return(splines_model)
}
