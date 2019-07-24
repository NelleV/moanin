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
create_splines_model = function(meta, formula=NULL, basis=NULL,
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
    } 

    splines_model = list()
    splines_model$degrees_of_freedom = degrees_of_freedom
    splines_model$"basis" = basis
    splines_model$"meta" = meta
    splines_model$"formula" = formula
    return(splines_model)
}
