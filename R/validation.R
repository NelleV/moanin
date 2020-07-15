



#' Internal Validation Checks
#' @keywords internal
#' @name internal
#' @return Does not return anything. Only hits errors if there are problems.
check_data_meta = function(data, object){
    if(!inherits(object,"Moanin")) stop("Internal coding error: passed object",
    "is not of class Moanin")
    dim_data = dim(data)
    dim_meta = dim(object)
    data = as.matrix(data)
    if(dim_meta[2] != dim_data[2]){
        stop(
            "User-given data and Moanin object are inconsistent. Data is has ", 
            ncol(data),"columns; Moanin object has", ncol(object))
    }

    if(!is.numeric(data)){
        stop("Data should be of type numeric")
    }

}


#' Check is 2D
#'
#' @keywords internal
#' @rdname internal
check_is_2d = function(X){
    dim_data = dim(X)
    if(is.null(dim_data)){
        stop("Data is expected to be 2D. No dimension found.")
    }
}


#' Check input are contrasts
#'
#' Will check that the contrasts provided are indeed contrasts. contrasts are
#'  expected to be either a vector of string or a matrix containing the
#'  contrasts coefficients.
#'  
#' @details If a vector of string is provided, the function will call
#'  limma::makeContrast in order to obtain the contrasts coefficients.
#'
#' @details If a contrasts matrix is provided, it will perform a number of
#'   checks on the contrasts matrix to make sure it contains the number of rows
#'   expected, and that each contrast indeed sums to 0.
#' @returns \code{is_contrasts} returns the contrasts, with any corrections.
#' @keywords internal
#' @rdname internal
is_contrasts = function(contrasts, moanin_model){
    if(!inherits(moanin_model, "Moanin") ) 
        stop("Coding error: internal function is_contrasts expect class",
            "Moanin object")
    if(is.vector(contrasts)){
        # XXX Should we add more tests here in order to provide meaningful
        # error messages?
        contrasts = limma::makeContrasts(
            contrasts=contrasts,
            levels=levels(group_variable(moanin_model)))
    }else{
        # Basic checks to make sure we have contrasts
        dim_contrasts = dim(contrasts)
        if(is.null(dim_contrasts)){
            msg = "Contrasts need to be either a vector of string or a 2D matrix"
            stop(msg)
        }
        
        n_groups = length(levels(group_variable(moanin_model)))
        if(dim(contrasts)[1] != n_groups){
            msg = paste(
                "When contrasts provided is a matrix, it should contain the",
                "same number of rows as the number of groups. Provided:",
                dim(contrasts[1]), "expected:", n_groups)
            stop(msg)
        }
        
        if(any(colSums(contrasts) != 0)){
            msg = paste("When contrasts provided is",
                "matrix, all columns should sum to 0.")
            stop(msg)
        }
    }
    return(contrasts) 
}
