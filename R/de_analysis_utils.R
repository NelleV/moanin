library("limma")

#' Estimates log fold change
#'
#' @param data The data in a matrix
#' @param moanin_model a moanin_model object, see \code{\link{create_moanin_model}}
#' @param contrasts The contrasts to consider
#' @param method method for calculating the log-fold change. See details. 
#' @details The following methods exist for calculating the log-fold change
#'   between conditions over time (default is "timecourse"):
#' \itemize{
#' \item{\code{timecourse}}{}
#' \item{\code{sum}}{}
#' \item{\code{max}}{}
#' \item{\code{timely}}{}
#' \item{\code{epicon}}{}
#' \item{\code{abs_sum}}{}
#' \item{\code{abs_squared_sum}}{}
#' \item{\code{min}}{}
#' }
#' @export
estimate_log_fold_change = function(data, moanin_model, 
				    contrasts, method=c("timecourse", "sum", "max", "timely", 
		    "epicon", "abs_sum", "abs_squared_sum", "min")){
    # Should check that data and meta is sorted identically
    meta = moanin_model$meta
    # Should check that the method is a known method
    method<-match.arg(method)
    contrasts = is_contrasts(contrasts, meta)

    if(method == "timely"){
	    log_fold_changes = lfc_per_time(data, meta, contrasts)
    }else if(method %in% c("sum", "max", "min", "abs_sum", "abs_squared_sum",
			   "timecourse", "epicon")){
	timely_lfc = lfc_per_time(data, meta, contrasts)
	timely_lfc_meta = reconstruct_meta_from_lfc(timely_lfc, split_char=":")
    	log_fold_changes = data.frame(row.names=row.names(data))
	for(contrast in colnames(contrasts)){
	    mask = (timely_lfc_meta$Group == contrast) & !is.na(colSums(timely_lfc))
	    if(method == "max"){
		log_fold_changes[, contrast] = rowMax(abs(timely_lfc[, mask]))
	    } else if(method == "min") {
		log_fold_changes[, contrast] = rowMin(abs(timely_lfc[, mask])) 
	    }else if(method == "abs_sum"){
		log_fold_changes[, contrast] = rowSums(abs(timely_lfc[, mask]))
	    }else if(method == "abs_squared_sum"){
		log_fold_changes[, contrast] = rowSums(timely_lfc[, mask]**2)
	    }else if(method == "epicon" || method == "timecourse"){
		log_fold_changes[, contrast] = (
		    rowMeans(abs(timely_lfc[, mask])) * sign(rowSums(timely_lfc[, mask])))
	    }else if(method == "sum"){
		log_fold_changes[, contrast] = rowSums(timely_lfc[, mask])
	    }
	}

    }

    return(log_fold_changes)
    
}


estimate_log_fold_change_sum = function(data, meta, contrasts){
    sample_coefficients = lapply(meta$Group, function(x) return(contrasts[x, ]))
    sample_coefficients = as.matrix(sample_coefficients)

    # First, do weekly contrasts

    row.names(sample_coefficients) = row.names(meta)
    log_fold_changes = data.frame(row.names=row.names(data))
    for(column in 1:ncol(sample_coefficients)){
	sample_coefficient = as.vector(unlist(sample_coefficients[, column]))
	log_fold_changes[, column] = as.matrix(data) %*% sample_coefficient 	
    }

}


data_summarize_per_time = function(data, meta){
    all_group_times = levels(meta$WeeklyGroup)

    log_fold_changes = data.frame(row.names=row.names(data))
    for(column in all_group_times){
	mask = meta$WeeklyGroup == column
	average_expr = rowSums(t(t(data) * mask))
	log_fold_changes[column] = average_expr 
    }
}


# XXX helper function to reconstruct metadat from adat
reconstruct_meta_from_lfc = function(data_per_time, split_char="."){
    meta_per_time = t(
	as.data.frame(strsplit(colnames(data_per_time), split_char, fixed=TRUE)))
    row.names(meta_per_time) = colnames(data_per_time)
    colnames(meta_per_time) = c("Group", "Timepoint")
    meta_per_time = as.data.frame(meta_per_time)
    meta_per_time[, "Timepoint"] = as.numeric(meta_per_time[, "Timepoint"])
    return(meta_per_time) 
}


lfc_per_time = function(data, meta, contrasts){
    meta$Timepoint = as.factor(meta$Timepoint)

    averaged_data = average_replicates(data, meta)
    averaged_meta = reconstruct_meta_from_lfc(averaged_data, split_char=":")

    averaged_meta$Timepoint = as.factor(averaged_meta$Timepoint)
    sample_coefficients = sapply(averaged_meta$Group, function(x) return(contrasts[x, ]))
    if(is.null(dim(sample_coefficients))){
	sample_coefficients = as.matrix(sample_coefficients)
    }else{
	sample_coefficients = t(sample_coefficients)
    }
    colnames(sample_coefficients) = colnames(contrasts)

    log_fold_changes = matrix(NA, dim(data)[1],
			      dim(contrasts)[2]*length(unique(meta$Timepoint)))
    row.names(log_fold_changes) = row.names(data)
    colnames(log_fold_changes) = sapply(
	colnames(sample_coefficients),
	function(x){sapply(unique(averaged_meta$Timepoint), function(t){paste0(x, ":", t)})})

    for(column in colnames(sample_coefficients)){
	
	sample_coefficient = as.vector(unlist(sample_coefficients[, column]))
	coef_data = t(t(averaged_data) * sample_coefficient)
	for(timepoint in unique(averaged_meta$Timepoint)){
	    mask = averaged_meta$Timepoint == timepoint
	    if(sum(mask) != dim(contrasts)[1]){
		next
	    }
	    colname = paste0(column, ":", as.character(timepoint))
	    log_fold_changes[, colname] = rowSums(coef_data[, mask])
	}
    }
    return(log_fold_changes)
}


average_replicates = function(data, meta){
    timepoint_group = droplevels(meta$Group:meta$Timepoint)
    all_levels = levels(timepoint_group)

    # selecting certain columns sometimes returns vectors and sometimes 
    # matrices depending on whether the user wants a single column or several
    # columns. rowMeans requires a matrix. Thus, here's a small implementation
    # that will work in all cases.
    .row_means = function(data){
	if(is.vector(data)){
	    return(data)
	}else{
	    return(rowMeans(data))
	}
    }

    replicate_averaged = sapply(all_levels,
				function(m){.row_means(data[, timepoint_group==m])})
    return(replicate_averaged)
}

    
