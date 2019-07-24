library("limma")
library("stats")
library("edgeR")


#' Fit weekly differential expression analysis
#'
#' @param data Gene expression data
#' @param splines_model splines_model
#' @param contrasts Contrast to use.
#' @param use_voom_weights boolean: whether to use voom weights or not
#'
#' @export
DE_timepoints = function(data, splines_model,
			 contrasts,
			 use_voom_weights=TRUE){
    meta = splines_model$meta

    design = stats::model.matrix(~WeeklyGroup + 0, data=meta)

    cleaned_colnames = gsub("WeeklyGroup", "", colnames(design))
    colnames(design) = cleaned_colnames

    allcontrasts = limma::makeContrasts(
        contrasts=contrasts,
        levels=design)

    if(use_voom_weights){
        y = edgeR::DGEList(counts=data)
	y = edgeR::calcNormFactors(y, method="upperquartile")
        v = limma::voom(y, design, plot=FALSE)
	v = limma::lmFit(v)
    }else{
	v = limma::lmFit(data, design)	
    }

    fit = limma::contrasts.fit(v, allcontrasts)
    fit = limma::eBayes(fit)
    contrast_names = colnames(fit$p.value)
    fit$adj.p.value = stats::p.adjust(fit$p.value, method="BH")
    dim(fit$adj.p.value) = dim(fit$p.value)
    colnames(fit$adj.p.value) = contrast_names

    combine_results = function(ii, fit2){
	contrast_formula = contrasts[ii]
	de_analysis = data.frame(row.names=row.names(data))

        base_colname = gsub(" ", "", contrast_formula, fixed=TRUE)
	colname_pval = paste(base_colname, "-pval", sep="")
	colname_qval = paste(base_colname, "-qval", sep="")
	colname_lfc = paste(base_colname, "-lfc", sep="")

        tt = limma::topTable(
            fit2, coef=ii, number=length(rownames(fit2$coef)),
            p.value=1, adjust.method="none",
            genelist=rownames(fit2$coef))
	de_analysis[colname_pval] = fit2$p.value[, contrast_formula]
	de_analysis[colname_qval] = fit2$adj.p.value[,  contrast_formula]
	de_analysis[colname_lfc] = tt$logFC
	return(de_analysis)
    }

    all_results = do.call("cbind",
			  lapply(1:length(contrast_names),
				 combine_results, fit2=fit))
    return(all_results)
}


#' Creates pairwise contrasts for all timepoints
#'
#' @param group1 First group to consider
#' @param group2 Second group to consider
#' @param splines_model splines model object
#'
#' @examples
#' data(shoemaker2015)
#' meta = shoemaker2015$meta
#' splines_model = moanin::create_splines_model(meta)
#' contrasts = create_timepoints_contrasts("C", "M", splines_model)
#'
#' @export
create_timepoints_contrasts = function(group1, group2, splines_model){
    meta = splines_model$meta
    meta = meta[meta$Group %in% c(group1, group2),]
    all_timepoints = sort(unique(meta$Timepoint))
    contrasts = rep(NA, length(all_timepoints))
    for(i in 1:length(all_timepoints)){
	# First, check that the two conditions have been sampled for this
	# timepoint
	timepoint = all_timepoints[i]
	submeta = meta[meta$Timepoint == timepoint, ]
	if(length(unique(submeta$WeeklyGroup)) == 2){
	    groups = as.character(unique(submeta$WeeklyGroup))
	    contrasts[i] = paste0(group1, ".", timepoint, "-", group2, ".", timepoint)
	}else if(length(unique(submeta$WeeklyGroup)) == 1){
	    if(unique(submeta$Group)[1] ==  group1){
		missing_condition = group2
	    }else{
		missing_condition = group1
	    }
	    msg = paste(
		"moanin::create_timepoints_contrasts: timepoint",
		timepoint, "is missing in condition", missing_condition)
	    warning(msg)
	}
    }
    return(contrasts[!is.na(contrasts)])
}
