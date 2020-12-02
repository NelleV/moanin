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
#' # compare groups within each timepoint
#' contrasts <- create_timepoints_contrasts(moanin,"C", "K",
#'    type="per_timepoint_group_diff")
#' head(contrasts)
#' deTimepoints=DE_timepoints(moanin, 
#'     contrasts=contrasts, use_voom_weights=FALSE)
#' head(deTimepoints)
#' # Control for replicate variable:
#' deTimepoints=DE_timepoints(moanin, 
#'     contrasts=contrasts, add_factors="Replicate",
#'     use_voom_weights=FALSE)
#' head(deTimepoints)
#'
#' # compare adjacent timepoints within each group
#' contrastsDiff <- create_timepoints_contrasts(moanin,"C",
#'    type="per_group_timepoint_diff")
#' deDiffTimepoints=DE_timepoints(moanin, 
#'     contrasts=contrastsDiff,
#'     use_voom_weights=FALSE)
#' # provide the sets of timepoints to compare:
#' contrastsDiff2<-create_timepoints_contrasts(moanin,"C",
#'    timepoints_before=c(72,120),timepoints_after=c(168,168),
#'    type="per_group_timepoint_diff")
#' deDiffTimepoints2=DE_timepoints(moanin, 
#'     contrasts=contrastsDiff2,
#'     use_voom_weights=FALSE)
#'
#' # Compare selected timepoints across groups. 
#' # This time we also return format="data.frame" which helps us keep track of
#' # the meaning of each contrast. 
#' contrastsGroupDiff<-create_timepoints_contrasts(moanin,"C", "K",
#'    timepoints_before=c(72,120),timepoints_after=c(168,168),
#'    type="group_and_timepoint_diff",format="data.frame")
#' head(contrastsGroupDiff)
#' deGroupDiffTimepoints=DE_timepoints(moanin, 
#'     contrasts=contrastsGroupDiff$contrasts,
#'     use_voom_weights=FALSE)
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
#'   that must match a value of the grouping variable contained in
#'   \code{moanin_model}.
#' @param group2 Second group to consider in making contrasts, character value
#'   that must match a value of the grouping variable contained in
#'   \code{moanin_model}, unless type=="per_group_timepoint_diff", in which case
#'   should be NULL (only \code{group1} is used in comparison)
#' @param timepoints vector of timepoints to compare. Must be contained in the
#'   \code{time_variable} of the \code{moanin} object.
#' @param timepoints_after for \code{type} equal to
#'   \code{"per_group_timepoint_diff"} or, \code{"group_and_timepoint_diff"},
#'   the set of timepoints to compare, see details. By default, taken from the
#'   \code{timepoints} variable.
#' @param timepoints_before for \code{type} equal to
#'   \code{"per_group_timepoint_diff"} or, \code{"group_and_timepoint_diff"},
#'   the set of timepoints to compare, see details. By default, taken from the
#'   \code{timepoints} variable.
#' @param format the choice of "vector" (the default) for
#'   \code{create_timepoints_contrasts} returns just the character vector of
#'   contrasts. If instead \code{format="data.frame"} then a data.frame is
#'   return that identifies the timepoint and group comparisons involved in each
#'   contrast. If this is the desired output, then the input to
#'   \code{DE_timepoints} should be the column corresponding to the contrast.
#'   See examples.
#' @param type the type of contrasts that should be created. See details.
#' @details \code{create_timepoints_contrasts} creates the needed contrasts for
#'   comparing groups or timepoints in the format needed for
#'   \code{DE_timepoints} (i.e. \code{\link[limma]{makeContrasts}}), to which the
#'   contrasts are ultimately passed. The time points and groups are determined
#'   by the levels of the \code{grouping_variable} and the values of
#'   \code{time_variable} in the \code{moanin_object} provided by the user.
#' @details Three different types of contrasts are created:
#'   \itemize{
#'  \item{"per_timepoint_group_diff"}{Contrasts that compare the groups within a
#'  timepoint}
#'  \item{"per_group_timepoint_diff"}{Contrasts that compare two timepoints
#'  within a group}
#'  \item{"group_and_timepoint_diff"}{Contrasts that compare the
#'   difference between two timepoints between two levels of the
#'   \code{group_variable} of the \code{Moanin} object. These are contrasts in
#'   the form (TP i - TP (i-1))[Group1] - (TP i - TP (i-1))[Group2].}
#' }
#' @export
#' @return \code{create_timepoints_contrasts}: a character vector with each
#'   element of the vector corresponding to a contrast to be compared.
#' @seealso \code{\link[limma]{makeContrasts}}
#' @rdname DE_timepoints
#' @importFrom utils head tail
#' @export
setMethod("create_timepoints_contrasts","Moanin",
 function(object, group1, group2=NULL, 
     type=c("per_timepoint_group_diff","per_group_timepoint_diff",
            "group_and_timepoint_diff"),
     timepoints=sort(unique(time_variable(object))),
     timepoints_before=head(sort(timepoints),-1),
     timepoints_after=tail(sort(timepoints),-1),
     format=c("vector","data.frame")
     ){
    type<-match.arg(type)
    format<-match.arg(format)
    if(type=="per_timepoint_group_diff"){
        if(is.null(group2)) 
            stop("cannot choose type='per_timepoint_group_diff'" + 
            "and give a NULL value for argument `group2`")
        if(!all(timepoints%in% time_variable(object))) 
            stop("timepoints must consist only of timepoints in the time_variable of Moanin object")
        
        contrasts<-pertimepoint_contrast(object=object, group1=group1,
            group2=group2,timepoints=timepoints)
    }
    if(type=="group_and_timepoint_diff"){
        if(is.null(group2)) 
            stop("cannot choose type='group_and_timepoint_diff'" + 
             "and give a NULL value for argument `group2`")
        if(!all(timepoints_before %in% time_variable(object))) 
            stop("timepoints_before must consist only of timepoints in the time_variable of Moanin object")
        if(!all(timepoints_after %in% time_variable(object))) 
            stop("timepoints_after must consist only of timepoints in the time_variable of Moanin object")
        contrasts<-timepointdiff_contrasts(object=object, group1=group1, 
            group2=group2, timepoints_before=timepoints_before,
            timepoints_after=timepoints_after)
    }
    if(type=="per_group_timepoint_diff"){
        if(!is.null(group2)) 
            stop("cannot choose type='per_group_timepoint_diff'" + 
                     "and give a value for argument `group2`")
        contrasts<-timepointdiff_contrasts(object=object, group1=group1, 
            group2=NULL, timepoints_before=timepoints_before,
            timepoints_after=timepoints_after)
    }
    if(format=="vector") return(contrasts$contrasts)
    else{
        return(contrasts)
    }
})
         
pertimepoint_contrast<-function(object, group1, group2, 
     timepoints){
    object <- object[,group_variable(object) %in% c(group1, group2)]
    all_timepoints <- sort(unique(time_variable(object)))
    timepoints<-timepoints[.which_timepoints(timepoints,all_timepoints,argname="timepoints")]
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
    timepoints<-timepoints[!is.na(contrasts)]
    contrasts<-contrasts[!is.na(contrasts)]
    return(data.frame("contrasts"=contrasts,"timepoints"=as.character(timepoints),"group"=paste0(group1,"-",group2)))
 }

.which_timepoints<-function(timepoints, possibles,argname){
    wh<-which(timepoints%in% possibles)
    if(!all(timepoints%in% possibles)){
        if(length(wh)>0){
            warning("removing timepoints in ",argname," not measured for these groups\n")
        }
        else{
            stop("None of the requested timepoints measured for these groups")
        }
    }
    return(wh)
}

timepointdiff_contrasts<-function(object, group1, group2, 
    timepoints_before=NULL,timepoints_after=NULL){
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
    wh_before<-.which_timepoints(timepoints_before,all_timepoints,"timepoints_before")
    wh_after<-.which_timepoints(timepoints_after,all_timepoints,"timepoints_after")
    wh<-intersect(wh_before,wh_after)
    timepoints_before<-timepoints_before[wh]
    timepoints_after<-timepoints_after[wh]
    
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
            msg <- paste(msg,"Cannot compare",tpbefore,"and",tpafter,
                "because one of the timepoints is missing in one of the conditions.\n")
            foundMissing<-TRUE
        }else{
          if(!is.null(group2))
              contrasts[i] <- paste0(group1, ".", tpafter, "-", group1, ".", tpbefore,"-",group2,".",tpafter,"+",group2,".",tpbefore)
          else
              contrasts[i] <- paste0(group1, ".", tpafter, "-", group1, ".", tpbefore)
        } 
    }
    if(foundMissing) warning(msg)
    timepoints_before<-timepoints_before[!is.na(contrasts)]
    timepoints_after<-timepoints_after[!is.na(contrasts)]
    if(!is.null(group2)) group<-paste0(group1,"-",group2) else group<-group1
    contrasts<-contrasts[!is.na(contrasts)]
    timepoints<-paste0(timepoints_after,"-",timepoints_before)
    return(data.frame("contrasts"=contrasts,"timepoints"=as.character(timepoints),"group"=group))
    
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
#' contrasts <- create_timepoints_contrasts(moanin, "C", "K")
#' deTimepoints <- DE_timepoints(moanin, 
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
            FUN=function(x){.subset2(x, 3)},
            character(1))
    }
    number_de_genes_per_time <- matrixStats::colSums2(
        de_results[, qval_colnames] < threshold)
    
    graphics::barplot(number_de_genes_per_time, 
            names.arg=labels, xlab=xlab,ylab=ylab,
            main=main, ...)  
}
