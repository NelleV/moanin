

#' Check that the metadata provided is what we expect
#'
#' This method will raise errors if the metadata provided is not as expected.
#'
#' @param meta metadata
#' @param check_replicates boolean, optional, default: FALSE
#'	  If TRUE, checks whether metadata contains a column Replicate
#' @return meta returns the metadata with additional columns if necessary.
#'
#' @keywords internal
check_meta = function(meta, check_replicates=FALSE){
    metadata_column_names = colnames(meta)
    if(!("Group" %in% metadata_column_names)){
	stop(
	    "moanin::create_moanin_model: " +
	    "Metadata doesn't contain expected information." +
	    " Group column is missing.")
    }

    if(!("Timepoint" %in% metadata_column_names)){
	stop(
	    "moanin::create_moanin_model: " +
	    "Metadata doesn't contain expected information." +
	    " Timepoint column is missing.")
    }

    if(check_replicates & !("Replicate" %in% metadata_column_names)){
	stop(
	    "moanin::create_moanin_model: " +
	    "Metadata doesn't contain expected information." +
	    "Replicate column is missing")
    }

    # Check that Timepoint is numeric.
    if(!is.numeric(meta$Timepoint)){
	stop(
	    "moanin::create_moanin_model: " +
	    "Timepoint column is expected to be numeric")
    }

    if(!is.factor(meta$Group)){
	stop(
	    "moanin::create_moanin_model: " +
	    "Group column is expected to be factors")
    }


    # Just create this one.
    if(!("WeeklyGroup" %in% metadata_column_names)){
	meta["WeeklyGroup"] = as.factor(
	    make.names(meta$Group:as.factor(meta$Timepoint)))
    }
    return(meta) 
}


#' Check data and meta
#' @keywords internal
check_data_meta = function(data, meta){
    dim_data = dim(data)
    dim_meta = dim(meta)
    data = as.matrix(data)
    if(dim_meta[1] != dim_data[2]){
	stop(
	    "Data and metadata are inconsistent. Data is of shape (Xx"+
	    "Metadata is of shape XX")
    }

    if(!is.numeric(data)){
	stop("Data should be of type numeric")
    }
    meta = check_meta(meta)
}


#' Check is 2D
#'
#' @keywords internal
check_is_2d = function(X){
    dim_data = dim(X)
    if(is.null(dim_data)){
	stop("Data is expected to be 2D. No dimension found.")
    }
}


#' Check is contrasts
#'
#' Will check that the contrasts provided are indeed contrasts. contrasts are
#'  expected to be either a vector of string or a matrix containing the
#'  contrasts coefficients.
#'  
#' If a vector of string is provided, the function will call
#'  limma::makeContrast in order to obtain the contrasts coefficients.
#'
#' If a contrasts matrix is provided, it will perform a number of checks on
#' the contrasts matrix to make sure it contains the number of rows expected,
#' and that each contrast indeed sums to 0.
#'
#' @keywords internal
is_contrasts = function(contrasts, meta){
    if(is.vector(contrasts)){
	# XXX Should we add more tests here in order to provide meaningful
	# error messages?
	contrasts = limma::makeContrasts(
	    contrasts=contrasts,
	    levels=levels(meta$Group))
    }else{
	# Basic checks to make sure we have contrasts
	dim_contrasts = dim(contrasts)
	if(is.null(dim_contrasts)){
	    msg = "Contrasts need to be either a vector of string or a 2D matrix"
	    stop(msg)
	}

	n_groups = length(levels(meta$Group))
	if(dim(contrasts)[1] != n_groups){
	    msg = paste(
		"When contrasts provided is a matrix, it should contain the",
		"same number of rows as the number of groups. Provided:",
		dim(contrats[1]), "expected:", n_groups)
	    stop(msg)
	}

	if(any(colSums(contrasts) != 0)){
	    msg = paste(
		"When contrasts provided is matrix, all columns should sum to 0.")
	    stop(msg)
	}
    }
    return(contrasts) 
}
