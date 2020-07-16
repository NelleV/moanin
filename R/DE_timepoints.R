setGeneric("DE_timepoints", 
           function(object,...) { standardGeneric("DE_timepoints")})
setGeneric("create_timepoints_contrasts", 
           function(object,...) { standardGeneric("create_timepoints_contrasts")})

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
#' @aliases create_timepoints_contrasts DE_timepoints,Moanin-method
#' @aliases create_timepoints_contrasts,Moanin-method
#' @name DE_timepoints
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom limma voom lmFit contrasts.fit eBayes
#' @details If \code{use_voom_weights=TRUE}, then the voom weights will be fit
#'   and the limma
#' @details If \code{use_voom_weights=TRUE}, the data is given directly to limma
#'   via \code{assay(object)}, unless log_transform(object)=TRUE, in which case
#'   it will be applied to \code{log(assay(object)+1)}.
#' @examples 
#' data(exampleData)
#' moanin = create_moanin_model(data=testData, meta=testMeta)
#' contrasts = create_timepoints_contrasts(moanin,"C", "K")
#' head(contrasts)
#' deTimepoints=DE_timepoints(moanin, 
#'     contrasts=contrasts, use_voom_weights=FALSE)
#' head(deTimepoints)
#' @export
setMethod("DE_timepoints","Moanin",
           function(object,
                         contrasts,
                         use_voom_weights=TRUE){
    
    design = stats::model.matrix(~WeeklyGroup + 0, data=colData(object))
    
    cleaned_colnames = gsub("WeeklyGroup", "", colnames(design))
    colnames(design) = cleaned_colnames
    
    allcontrasts = limma::makeContrasts(
        contrasts=contrasts,
        levels=design)
    
    if(use_voom_weights){
        y = edgeR::DGEList(counts=assay(object))
        y = edgeR::calcNormFactors(y, method="upperquartile")
        v = limma::voom(y, design, plot=FALSE)
        v = limma::lmFit(v)
    }else{
        y<-get_log_data(object)
        v = limma::lmFit(y, design)
    }
    
    fit = limma::contrasts.fit(v, allcontrasts)
    fit = limma::eBayes(fit)
    contrast_names = colnames(fit$p.value)
    fit$adj.p.value = stats::p.adjust(fit$p.value, method="BH")
    dim(fit$adj.p.value) = dim(fit$p.value)
    colnames(fit$adj.p.value) = contrast_names
    
    combine_results = function(ii, fit2){
        contrast_formula = contrasts[ii]
        de_analysis = data.frame(row.names=row.names(object))
        
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
)

#' Creates pairwise contrasts for all timepoints
#'
#' @param group1 First group to consider in making contrasts, character value
#'   that must match a value contained in \code{moanin_model$meta}.
#' @param group2 Second group to consider in making contrasts, character value
#'   that must match a value contained in \code{moanin_model$meta}.
#' @details \code{create_timepoints_contrasts} creates the needed contrasts for
#'   comparing two groups for every timepoint in the format needed for
#'   \code{DE_timepoints} (i.e. \code{\link[limma]{makeContrasts}}, to which the
#'   contrasts are ultimately passed). The time points are determined by the
#'   meta data in the \code{moanin_object} provided by the user.
#' @return \code{create_timepoints_contrasts}: a character vector with each
#'   element of the vector corresponding to a contrast to be compared.
#' @seealso \code{\link[limma]{makeContrasts}}
#' @rdname DE_timepoints
#' @export
setMethod("create_timepoints_contrasts","Moanin",
 function(object, group1, group2){
    object = object[,group_variable(object) %in% c(group1, group2)]
    all_timepoints = sort(unique(time_variable(object)))
    contrasts = rep(NA, length(all_timepoints))
    msg<-""
    foundMissing<-FALSE
    for(i in seq_along(all_timepoints)){
        # First, check that the two conditions have been sampled for this
        # timepoint
        timepoint = all_timepoints[i]
        submeta = object[,time_variable(object) == timepoint]
        if(length(unique(time_by_group_variable(submeta))) == 2){
            groups = as.character(unique(time_by_group_variable(submeta)))
            contrasts[i] = paste0(group1, ".", timepoint, "-", group2, ".", 
                                  timepoint)
        }else if(length(unique(time_by_group_variable(submeta))) == 1){
            if(unique(group_variable(submeta))[1] ==  group1){
                missing_condition = group2
            }else{
                missing_condition = group1
            }
            msg = paste0(msg,paste("timepoint",
                                   timepoint, "is missing in condition", 
                                   missing_condition,"\n"))
            foundMissing<-TRUE
        }
    }
    if(foundMissing) warning(msg)
    return(contrasts[!is.na(contrasts)])
 }
)

#' Creates pairwise contrasts for all timepoints
#'
#' @param de_results results from \code{\link{DE_timepoints}}
#' @param type type of p-value to count ("qval" or "pval")
#' @param labels labels to give each bar
#' @param threshold cutoff for counting gene as DE
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param main title of plot
#' @param ... arguments passed to \code{\link{barplot}}
#' @details \code{create_timepoints_contrasts} creates the needed contrasts for
#'   comparing two groups for every timepoint in the format needed for
#'   \code{DE_timepoints} (i.e. \code{\link[limma]{makeContrasts}}, to which the
#'   contrasts are ultimately passed). The time points are determined by the
#'   meta data in the \code{moanin_object} provided by the user.
#' @return This is a plotting function, and returns (invisibly) the results of 
#'   \code{\link{barplot}}
#' @aliases perWeek_barplot
#' @examples 
#' data(exampleData)
#' moanin = create_moanin_model(data=testData, meta=testMeta)
#' contrasts = create_timepoints_contrasts(moanin,"C", "K")
#' deTimepoints=DE_timepoints(moanin, 
#'     contrasts=contrasts, use_voom_weights=FALSE)
#' perWeek_barplot(deTimepoints)
#' @export
perWeek_barplot = function(de_results, type=c("qval","pval"),
                          labels=NULL, threshold=0.05,
                          xlab="Timepoint", ylab="Number of DE genes", main="", 
                           ...){
    
    type<-match.arg(type)
    qval_colnames = colnames(de_results)[
        grepl(type, colnames(de_results))]
    if(is.null(labels)){
        stringReplace<-paste0("_",type)
        labels = sapply(
            strsplit(gsub(stringReplace, "", qval_colnames), "\\."), .subset2, 3)
    }
    number_de_genes_per_time = colSums(de_results[, qval_colnames] < threshold)
    
    graphics::barplot(number_de_genes_per_time, 
            names.arg=labels, xlab=xlab,ylab=ylab,
            main=main, ...)  
}
