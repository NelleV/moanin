# Creates a "class" containing the information related to the splines used.

#' Create a splines model
#'
#' This model is used to conserve information relating to the model, such as
#' the formula, the basis, the meta data.
#'
#' @param meta	data.frame containing the metadata. Metadata needs to contain
#'	column "Group" and "Timepoint"
#' @param formula formula object, optional, default: NUlL
#'	Formula of the splines, e.g:
#'	    formula = ~Group:splines::ns(Timepoint, df=4) + Group + 0
#' @param basis	matrix, optional, default: NULL
#'	A basis object.
#' @param degrees_of_freedom int, optional, default: 4
#'	Number of degrees of freedom to use if neither the basis nor the
#'	formula is provided	 
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
#' @param meta	data.frame containing the metadata. Metadata needs to contain
#'	column "Group" and "Timepoint"
#' @param formula formula object, optional, default: NUlL
#'	Formula of the splines, e.g:
#'	    formula = ~Group:splines::ns(Timepoint, df=4) + Group + 0
#' @param basis	matrix, optional, default: NULL
#'	A basis object.
#' @param degrees_of_freedom int, optional, default: 4
#'	Number of degrees of freedom to use if neither the basis nor the
#'	formula is provided	 
#'
#' @return An object containing the following information:
#'
#'	- The basis: a matrix of shape (degrees_of_freedom, n_samples)
#'	- The metadata: a data.frame of shape (n_samples, n_metadat)
#'	- The formula used, when provided.
#'	- The number of degrees of freedom.
#'
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
