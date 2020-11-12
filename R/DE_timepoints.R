setGeneric("DE_timepoints", 
           function(object,...) { standardGeneric("DE_timepoints")})
setGeneric("create_timepoints_contrasts", 
           function(object,...) { standardGeneric("create_timepoints_contrasts")})
setGeneric("create_diff_contrasts", 
                      function(object,...) { standardGeneric("create_diff_contrasts")})


#' Fit weekly differential expression analysis
#'
#' @inheritParams DE_timecourse
#' @param add_factors A character vector of additional variables to add to the 
#' design. See details. 
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
#' @details By default the formula fitted for each gene is
#' \preformatted{   
#'    ~ Group*Timepoint +0
#' }
#' If the user gives values to \code{add_factors}, then the vector of character 
#' values given in \code{add_factors} will be \emph{added} to the default formula. 
#' So that \code{add_factors="Replicate"} will change the formula to
#' \preformatted{   
#'    ~ Group*Timepoint +0 + Replicate
#' }
#' This allows for a small amount of additional complexity to control 
#' for other variables. Users should work directly with limma for 
#' more complex models. 
#' @details If \code{use_voom_weights=TRUE}, the data is given directly to limma
#'   via \code{assay(object)}. The specific series of
#'   calls is:
#' \preformatted{   
#'    y <- edgeR::DGEList(counts=assay(object))
#'    y <- edgeR::calcNormFactors(y, method="upperquartile")
#'    v <- limma::voom(y, design, plot=FALSE)
#'    v <- limma::lmFit(v) 
#'    }
#' @details If the user set \code{log_transform=TRUE} in the creation of the
#'   \code{Moanin} object, this will not have an impact in the analysis if
#'   \code{use_voom_weights=TRUE}. Only if \code{use_voom_weights=FALSE} will
#'   this matter, in which case the log of the input data will be given to a
#'   regular call to \code{limma}:
#' \preformatted{
#'    y<-get_log_data(object)
#'    v <- limma::lmFit(y, design)
#' }
#' @examples 
#' data(exampleData)
#' moanin <- create_moanin_model(data=testData, meta=testMeta)
#' contrasts <- create_timepoints_contrasts(moanin,"C", "K")
#' head(contrasts)
#' deTimepoints=DE_timepoints(moanin, 
#'     contrasts=contrasts, use_voom_weights=FALSE)
#' head(deTimepoints)
#' contrastsDiff <- create_diff_contrasts(moanin,"C", "K")
#' deDiffTimepoints=DE_timepoints(moanin, 
#'     contrasts=contrasts,
#'     use_voom_weights=FALSE)
#' # provide the sets of timepoints to compare:
#' contrastsDiff<-create_timepoints_contrasts(moanin,"C", "K",timepoints_before=c(72,120),timepoints_after=c(168,168),type="diff_timepoint")
#' # Control for replicate variable:
#' deTimepoints=DE_timepoints(moanin, 
#'     contrasts=contrasts, add_factors="Replicate",
#'     use_voom_weights=FALSE)
#' head(deTimepoints)

#' @export
setMethod("DE_timepoints","Moanin",
           function(object,
                    contrasts,add_factors=NULL,
                    use_voom_weights=TRUE){
    designText<-"~WeeklyGroup + 0"
    if(!is.null(add_factors)){
        designText<-paste(designText,"+",
            paste(add_factors,collapse="+"))
    }
    designFormula<-stats::as.formula(designText)
    design <- try(stats::model.matrix(designFormula, data=colData(object)), silent=TRUE)
    if( inherits(design, "try-error")){
        stop("Error in creating the design matrix. Error:\n",design)
    }
    
    cleaned_colnames <- gsub("WeeklyGroup", "", colnames(design))
    colnames(design) <- cleaned_colnames
    
    allcontrasts <- limma::makeContrasts(
        contrasts=contrasts,
        levels=design)
    
    if(use_voom_weights){
        y <- edgeR::DGEList(counts=assay(object))
        y <- edgeR::calcNormFactors(y, method="upperquartile")
        v <- limma::voom(y, design, plot=FALSE)
        v <- limma::lmFit(v)
    }else{
        y<-get_log_data(object)
        v <- limma::lmFit(y, design)
    }
    
    fit <- limma::contrasts.fit(v, allcontrasts)
    fit <- limma::eBayes(fit)
    contrast_names <- colnames(fit$p.value)
    fit$adj.p.value <- stats::p.adjust(fit$p.value, method="BH")
    dim(fit$adj.p.value) <- dim(fit$p.value)
    colnames(fit$adj.p.value) <- contrast_names
    
    combine_results <- function(ii, fit2){
        contrast_formula <- contrasts[ii]
        de_analysis <- data.frame(row.names=row.names(object))
        
        base_colname <- gsub(" ", "", contrast_formula, fixed=TRUE)
        colname_pval <- paste(base_colname, "_pval", sep="")
        colname_qval <- paste(base_colname, "_qval", sep="")
        colname_lfc <- paste(base_colname, "_lfc", sep="")
        
        tt <- limma::topTable(
            fit2, coef=ii, number=length(rownames(fit2$coef)),
            p.value=1, adjust.method="none",
            genelist=rownames(fit2$coef))
        de_analysis[colname_pval] <- fit2$p.value[, contrast_formula]
        de_analysis[colname_qval] <- fit2$adj.p.value[,  contrast_formula]
        de_analysis[colname_lfc] <- tt$logFC
        return(de_analysis)
    }
    
    all_results <- do.call("cbind",
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
#' @param timepoints_after for create_diff_contrasts, the set of timepoints to compare, see details where contrast will be timepoints_after - timepoints_before 
#' @param timepoints_before for create_diff_contrasts, see details the set of timepoints to compare, where contrast will be timepoints_after - timepoints_before 
#' @details \code{create_timepoints_contrasts} creates the needed contrasts for
#'   comparing two groups for every timepoint in the format needed for
#'   \code{DE_timepoints} (i.e. \code{\link[limma]{makeContrasts}}, to which the
#'   contrasts are ultimately passed). The time points are determined by the
#'   meta data in the \code{moanin_object} provided by the user.
#' @details \code{type="diff_timepoint"} will create contrasts that compare the difference between two timepoints between two levels of the \code{group_variable} of the \code{Moanin} object. These are contrasts in the form (TP i - TP (i-1))[Group1] - (TP i - TP (i-1))[Group2]. 
#' @export
#' @return \code{create_timepoints_contrasts}: a character vector with each
#'   element of the vector corresponding to a contrast to be compared.
#' @seealso \code{\link[limma]{makeContrasts}}
#' @rdname DE_timepoints
#' @export
setMethod("create_timepoints_contrasts","Moanin",
 function(object, group1, group2, 
     type=c("per_timepoint","diff_timepoint"),
     timepoints=sort(unique(time_variable(object))),
     timepoints_before=head(sort(timepoints),-1),
     timepoints_after=tail(sort(timepoints),-1)
     ){
    type<-match.arg(type)
    if(type=="per_timepoint"){
        return(pertimepoint_contrast(object=object, group1=group1,
            group2=group2,timepoints=timepoints))
    }
    if(type=="diff_timepoint"){
        return(timepointdiff_contrasts(object=object, group1=group1, 
            group2=group2, timepoints_before=timepoints_before,
            timepoints_after=timepoints_after))
    }
})
         
pertimepoint_contrast<-function(object, group1, group2, 
     timepoints){
    object <- object[,group_variable(object) %in% c(group1, group2)]
    all_timepoints <- sort(unique(time_variable(object)))
    if(!all(timepoints%in% all_timepoints)) 
        stop("timepoints must consist only of timepoints in the time_variable of Moanin object")
    contrasts <- rep(NA, length(timepoints))
    msg<-""
    foundMissing<-FALSE
    for(i in seq_along(timepoints)){
        # First, check that the two conditions have been sampled for this
        # timepoint
        timepoint <- timepoints[i]
        submeta <- object[,time_variable(object) == timepoint]
        if(length(unique(time_by_group_variable(submeta))) == 2){
            contrasts[i] <- paste0(group1, ".", timepoint, "-", group2, ".", 
                                  timepoint)
        }else if(length(unique(time_by_group_variable(submeta))) == 1){
            if(unique(group_variable(submeta))[1] ==  group1){
                missing_condition <- group2
            }else{
                missing_condition <- group1
            }
            msg <- paste0(msg,paste("timepoint",
                                   timepoint, "is missing in condition", 
                                   missing_condition,"\n"))
            foundMissing<-TRUE
        }
    }
    if(foundMissing) warning(msg)
    return(contrasts[!is.na(contrasts)])
 }


timepointdiff_contrasts<-function(object, group1, group2, timepoints_before=NULL,timepoints_after=NULL){
    object <- object[,group_variable(object) %in% c(group1, group2)]
    all_timepoints <- sort(unique(time_variable(object)))
    
    ### Checks for timepoints
    if((is.null(timepoints_before) & !is.null(timepoints_after)) ||
    (!is.null(timepoints_before) & is.null(timepoints_after))){
        stop("either timepoints_before and timepoints_after must be given, or both must be NULL")
    }
    if(is.null(timepoints_before)){
        timepoints_before<-head(all_timepoints,-1)
        timepoints_after<-tail(all_timepoints,-1)
    }
    if(!all(timepoints_before %in% all_timepoints) || !all(timepoints_after %in% all_timepoints)) 
        stop("timepoints_before and timepoints_after must consist only of timepoints in the time_variable of Moanin object")
    if(!all(timepoints_before<timepoints_after)) 
        stop("each timepoints_after element has to be strictly greater than the corresponding timepoints_before element")
    
    contrasts <- rep(NA, length(timepoints_before))
    # Will give a tally of timepoint pairs can't do
    msg<-""
    foundMissing<-FALSE
    for(i in seq_along(timepoints_before)){
        tpbefore <- timepoints_before[i]
        tpafter<-timepoints_after[i]
        
        # First, check that the two conditions have been sampled 
        # for both timepoints
        # Could do all at once, but not worth effort
        combos<-expand.grid(tp=c(tpbefore,tpafter),
            groups=c(group1,group2))
        npercombo<-sapply(1:nrow(combos),function(i){
            tp<-combos[i,1]
            gp<-as.character(combos[i,2])
            sum(time_variable(object)==tp & group_variable(object)==gp)
        })
        if(any(npercombo==0)){
            msg <- paste(msg,"Cannot compare",tpbefore,"and",tpafter,"because one of the timepoints is missing in one of the conditions.\n")
            foundMissing<-TRUE
        }else{
          contrasts[i] <- paste0(group1, ".", tpafter, "-", group1, ".", tpbefore,"-",group2,".",tpafter,"+",group2,".",tpbefore)
            
        } 
    }
    if(foundMissing) warning(msg)
    return(contrasts[!is.na(contrasts)])
 }


#' Creates barplot of results of per-timepoint comparison
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
#' moanin <- create_moanin_model(data=testData, meta=testMeta)
#' contrasts <- create_timepoints_contrasts(moanin,"C", "K")
#' deTimepoints=DE_timepoints(moanin, 
#'     contrasts=contrasts, use_voom_weights=FALSE)
#' perWeek_barplot(deTimepoints)
#' @export
perWeek_barplot <- function(de_results, type=c("qval","pval"),
                          labels=NULL, threshold=0.05,
                          xlab="Timepoint", ylab="Number of DE genes", main="", 
                           ...){
    
    type <- match.arg(type)
    qval_colnames <- colnames(de_results)[
        grepl(type, colnames(de_results))]
    if(is.null(labels)){
        stringReplace <- paste0("_",type)
        labels <- vapply(
            strsplit(gsub(stringReplace, "", qval_colnames), "\\."),
            .subset2, 3,
            character(1))
    }
    number_de_genes_per_time <- matrixStats::colSums2(
        de_results[, qval_colnames] < threshold)
    
    graphics::barplot(number_de_genes_per_time, 
            names.arg=labels, xlab=xlab,ylab=ylab,
            main=main, ...)  
}
