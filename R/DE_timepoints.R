
#' Fit weekly differential expression analysis
#'
#' @inheritParams DE_timecourse
#' @return A \code{data.frame} with three columns for each of the contrasts
#'   given in \code{contrasts}, corresponding to the raw p-value of the contrast
#'   for that gene (\code{_pval}), the adjusted p-value (\code{_qval}), and the
#'   estimate of log-fold-change (\code{_lfc}). The adjusted p-values are
#'   FDR-adjusted based on the Benjamini-Hochberg method, as implemented in
#'   \code{\link[stats]{p.adjust}}. The adjustment is done across all p-values
#'   for all contrasts calculated.
#' @aliases create_timepoints_contrasts
#' @examples 
#' data(exampleData)
#' moanin = create_moanin_model(testMeta)
#' contrasts = create_timepoints_contrasts("C", "K", moanin)
#' head(contrasts)
#' deTimepoints=DE_timepoints(data=testData, moanin_model=moanin, 
#'     contrasts=contrasts, use_voom_weights=FALSE)
#' head(deTimepoints)
#' @export
DE_timepoints = function(data, moanin_model,
                         contrasts,
                         use_voom_weights=TRUE){
    meta = moanin_model$meta
    
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
        colname_pval = paste(base_colname, "_pval", sep="")
        colname_qval = paste(base_colname, "_qval", sep="")
        colname_lfc = paste(base_colname, "_lfc", sep="")
        
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
                          lapply(seq_along(contrast_names),
                                 combine_results, fit2=fit))
    return(all_results)
}


#' Creates pairwise contrasts for all timepoints
#'
#' @param group1 First group to consider in making contrasts, character value that must match a
#'   value contained in \code{moanin_model$meta}.
#' @param group2 Second group to consider in making contrasts, character value that must match a
#'   value contained in \code{moanin_model$meta}.
#' @details \code{create_timepoints_contrasts} creates the needed contrasts for comparing two groups for every
#'   timepoint in the format needed for \code{DE_timepoints} (i.e.
#'   \code{\link[limma]{makeContrasts}}, to which the contrasts are ultimately
#'   passed). The time points are determined by the meta data in the
#'   \code{moanin_object} provided by the user.
#' @return \code{create_timepoints_contrasts}: a character vector with each element of the vector corresponding to a
#'   contrast to be compared.
#' @seealso \code{\link[limma]{makeContrasts}}
#' @rdname DE_timepoints
#' @export
create_timepoints_contrasts = function(group1, group2, moanin_model){
    meta = moanin_model$meta
    gpVar = moanin_model$group_variable
    tpVar = moanin_model$time_variable
    meta = meta[meta[,gpVar] %in% c(group1, group2),]
    all_timepoints = sort(unique(meta[,tpVar]))
    contrasts = rep(NA, length(all_timepoints))
    msg<-""
    foundMissing<-FALSE
    for(i in seq_along(all_timepoints)){
        # First, check that the two conditions have been sampled for this
        # timepoint
        timepoint = all_timepoints[i]
        submeta = meta[meta[,tpVar] == timepoint, ]
        if(length(unique(submeta$WeeklyGroup)) == 2){
            groups = as.character(unique(submeta$WeeklyGroup))
            contrasts[i] = paste0(group1, ".", timepoint, "-", group2, ".", timepoint)
        }else if(length(unique(submeta$WeeklyGroup)) == 1){
            if(unique(submeta[,gpVar])[1] ==  group1){
                missing_condition = group2
            }else{
                missing_condition = group1
            }
            msg = paste0(msg,paste("timepoint",
                                   timepoint, "is missing in condition", missing_condition,"\n"))
            foundMissing<-TRUE
        }
    }
    if(foundMissing) warning(msg)
    return(contrasts[!is.na(contrasts)])
}
