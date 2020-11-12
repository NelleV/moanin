setGeneric("DE_timecourse", 
           function(object,...) { standardGeneric("DE_timecourse")})

# This file contains EDGE code adapted to work with limma contrasts

center_data <- function(y, ng_labels){
    for(g in levels(ng_labels)){
        whKeep <- which(ng_labels == g)
        sub_mean <- matrixStats::rowMeans2(y[, whKeep])
        y[, whKeep] <- y[, whKeep] - sub_mean
    }
    return(y)
}

compute_beta_null <- function(basis, beta, contrasts_coef){
    ng <- length(contrasts_coef)
    df <- ncol(basis) / ng
    contrasts_coef_ <- rep(contrasts_coef, times=df)
    
    # Reshape b so that each row corresponds to a group and drop the intercept
    b_ <- array(beta, dim=c(dim(beta)[1], ng, df))
    # First start by constructing the matrix T
    t_ <- apply(b_, 1, function(x){matrixStats::colSums2(x*contrasts_coef)})
    
    # The observations can not be assumed to be balanced...
    # We need to get rid of the intercept for this part
    # FIXME don't invert this matrix…
    part_K <- MASS::ginv(t(basis) %*% basis)
    K <- part_K * contrasts_coef_**2
    
    # We now need to sum all elements associated to the same pairs of splines.
    # Which is, in a particular case every four elements in both directions.
    
    K <- vapply(seq_len(df),
                function(jg) matrixStats::rowSums2(K[, (jg-1)*ng + seq_len(ng)]),
                numeric(ncol(K)))
    K <- vapply(seq_len(df),
                function(jg) matrixStats::colSums2(K[(jg-1)*ng + seq_len(ng),]),
                numeric(ncol(K)))
    
    T_ <- MASS::ginv(K) %*% t_
    
    # We got T. Now, let's move on to the rest
    tmp <- as.array(rep(as.vector(T_), each=ng), dim=c(1, 1, 1))
    dim(tmp) <- c(ng * df, dim(beta)[1])
    C_ <- contrasts_coef * part_K %*% tmp
    
    C_ <- t(C_)
    dim(C_) <- c(dim(beta)[1], ng, df)
    beta_null <- b_ - C_
    
    # Last step: reshape beta null so that it is of the same shape as beta
    dim(beta_null) <- dim(beta)
    return(beta_null)
}


lrtStat <- function(resNull, resFull, ng_labels=NULL) {
    # FIbasisME I'm pretty sure that in the case of contrasts, the degrees of
    # freedom computed here are wrong as they include part of the data that is
    # not used for the test. This needs to be fixed
    stat <- 0
    if(is.null(ng_labels)){
        ss0 <- matrixStats::rowSums2(resNull^2)
        ss1 <- matrixStats::rowSums2(resFull^2)
        n <- ncol(resNull)
        stat <- stat + n * (ss0 - ss1)/(ss1)
        
    }else{
        for(g in levels(ng_labels)){
            whKeep <- which(ng_labels == g)
            sub_resNull <- resNull[,whKeep]
            sub_resFull <- resFull[,whKeep]
            
            # Somehow the two lines above don't return the same object depending on
            # the dimension of resNull and resFull, so need to distinguish the case
            # where there is only one observation in data.
            if(is.null(dim(sub_resNull))){
                ss0 <- sum(sub_resNull^2)
                ss1 <- sum(sub_resFull^2)
                n <- length(sub_resNull)
            }else{
                ss0 <- matrixStats::rowSums2(sub_resNull^2)
                ss1 <- matrixStats::rowSums2(sub_resFull^2)
                n <- ncol(sub_resNull)
            }
            stat <- stat + n * (ss0 - ss1)/(ss1)
        }
    }
    
    return(stat)
}

compute_pvalue <- function(basis, y, beta, beta_null, ng_labels,
                          n_groups=NULL,
                          n_samples=NULL,
                          degrees_of_freedom=NULL,
                          statistics="lrt",
                          df2=NULL, weights=NULL){
    if(inherits(y,"DataFrame")){
        y <- data.matrix(y)
    }
    
    fitFull <- beta %*% t(basis)
    
    fitNull <- beta_null %*% t(basis)
    
    if(!is.null(weights)){
        resNull <- weights^(1/2) * (y - fitNull)
        resFull <- weights^(1/2) * (y - fitFull)
    }else{
        resNull <- y - fitNull
        resFull <- y - fitFull
    }
    
    # estimate degrees of freedom.
    if(is.null(n_groups)){
        n_groups <- nlevels(ng_labels)
        # FIbasisME Raise warning
    }
    
    if(is.null(n_samples)){
        n_samples <- ncol(basis)
        # FIbasisME raise warning
    }
    
    df <- degrees_of_freedom
    
    if(statistics == "ftest"){
        stat <- lrtStat(resNull, resFull)
        if(is.null(df2)){
            # FIbasisME Check this.
            df2 <- n_samples - degrees_of_freedom * n_groups
        }
        df1 <- df
        pval <- stats::pf(stat * df2 / df1, df1=df1, df2=df2, lower.tail=FALSE)
    }else{
        lstat <- lrtStat(resNull, resFull, ng_labels=ng_labels)
        pval <- stats::pchisq(lstat, df=degrees_of_freedom, lower.tail=FALSE)
    }
    return(pval)
}

summarise <- function(basis, ng_levels) {
    basis_mean <- matrix(nrow=nrow(basis), ncol=nlevels(ng_levels))
    colnames(basis_mean) <- levels(ng_levels)
    rownames(basis_mean) <- rownames(basis)
    
    for(g in levels(ng_levels)){
        whKeep <- which(ng_levels == g)
        if(length(whKeep) > 1){
            basis_mean[, g] <- matrixStats::rowMeans2(basis[, whKeep])
        }else if(length(whKeep) != 0){
            basis_mean[, g] <- basis[, whKeep]
        }
    }
    return(basis_mean)
}


#' Run spline models and test for DE of contrasts.
#' 
#' @param object An object of class \code{\link{Moanin}}, an object containing
#'   all related information for time course data and the splines model that
#'   will be used (if applicable). See \code{\link{create_moanin_model}} for
#'   more details.
#'@param contrasts Contrasts, either provided as a vector of strings, or a
#'  matrix of contrasts coefficients obtained using
#'  \code{\link[limma]{makeContrasts}} from the package \code{limma}. If given
#'  as a character string, will be passed to \code{\link[limma]{makeContrasts}}
#'  to be converted into such a matrix.
#' @param center boolean, whether to center the data matrix
#' @param use_voom_weights boolean, optional, default: TRUE. 
#'  Whether to use voom weights. See details.
#' @details The implementation of the spline fit and the calculation of p-values
#'   was based on code from \code{\link[edge]{edge}}, and expanded to enable
#'   handling of comparisons of groups via contrasts. The code assumes that the \code{Moanin} object was created via either a formula or a basis where a different spline was fit for each \code{group_variable} and thus the contrasts are comparisons of those spline fits. If the \code{Moanin} object was created via user-provided basis matrix or formula, then the user should take a great deal of caution in using this code, as the degrees of freedom for the tests of significance cannot be verified to be correct. 
#' @seealso
#'   \code{\link[limma]{makeContrasts}}, \code{\link{create_moanin_model}},
#'   \code{\link{DE_timepoints}}, \code{\link[edge]{edge}}
#' @details If \code{use_voom_weights=TRUE}, then before fitting splines to each gene,
#' voom weights are calculated from \code{assay(object)}:
#' \preformatted{
#'    y <- edgeR::DGEList(counts=assay(object))
#'    y <- edgeR::calcNormFactors(y, method="upperquartile")
#'    v <- limma::voom(y, design, plot=FALSE)
#'    weights <- v$weights
#' }
#' The design matrix for the voom weights is based on the formula
#' \code{~Group + Timepoint +0}
#' where Group and Timepoint are replaced with the user-defined values where appropriate. 
#' These weights are given to the \code{lm.fit} which fits the spline coefficients.
#' This workflow assumes that the input to the \code{Moanin} object were counts.
#' @details If the user set \code{log_transform=TRUE} in the creation of the
#'   \code{Moanin} object, the splines will be fit to the log of the input data,
#'   and not directly to the input data. This is independent of whether the user
#'   chooses \code{use_voom_weights}.
#' @return A \code{data.frame} with two columns for each of the contrasts given
#'   in \code{contrasts}, corresponding to the raw p-value of the contrast for
#'   that gene (\code{_pval}) and the adjusted p-value (\code{_qval}). The
#'   adjusted p-values are FDR-adjusted based on the Benjamini-Hochberg method,
#'   as implemented in \code{\link[stats]{p.adjust}}. The adjustment is done
#'   across all p-values for all contrasts calculated.
#' @examples 
#' data(exampleData)
#' moanin <- create_moanin_model(data=testData, meta=testMeta)
#' deTimecourse=DE_timecourse(moanin, 
#'    contrasts="K-C", use_voom_weights=FALSE)
#' head(deTimecourse)
#' @name DE_timecourse
#' @aliases DE_timecourse,Moanin-method
#' @export
setMethod("DE_timecourse","Moanin",
         function(object,
                  contrasts,
                  center=FALSE,
                  use_voom_weights=TRUE){
    basis <- basis_matrix(object)
    
    ng_labels <- group_variable(object)
    ng <- nlevels(ng_labels)

    contrasts <- is_contrasts(contrasts, object)
    
    if(use_voom_weights){
        formulaText<-paste0("~",group_variable_name(object)," + ",
                            time_variable_name(object)," + 0")
        voomFormula = stats::as.formula(formulaText)
        
        design <- stats::model.matrix(voomFormula,
            data=colData(object))
        
        y <- edgeR::DGEList(counts=assay(object))
        y <- edgeR::calcNormFactors(y, method="upperquartile")
        v <- limma::voom(y, design, plot=FALSE)
        weights <- v$weights
    }else{
        weights <- NULL
    }
    
    y<-get_log_data(object)
    
    if(center){
        y <- center_data(y)
        basis <- t(center_data(t(basis)))
    }
    
    beta <- fit_splines(data=y, moanin_model=object, weights=weights)
    
    if(dim(contrasts)[1] != ng){
        stop("The contrast coef vector should be of the same size" +
                 " as the number of groups")
    }
    results <- data.frame(row.names=row.names(object))
    for(col in seq_len(ncol(contrasts))){
        contrast <- contrasts[, col]
        
        # Create the name of the column
        contrast_name <- colnames(contrasts)[col]
        contrast_name <- gsub(" ", "", contrast_name, fixed=TRUE)
        
        # Get the number of samples used for this particular contrast:
        groups_of_interest <- names(contrast)[contrast != 0]
        n_samples_fit <- sum(group_variable(object) %in% groups_of_interest)
        n_groups <- length(groups_of_interest)
        
        beta_null <- compute_beta_null(basis, beta, contrast)
        
        ## FIXME: This assumes a particular form for the formula, which may not be true if user add additional controlling values, for example. 
        degrees_of_freedom <- dim(basis)[2] / ng    
        pval <- compute_pvalue(basis, y, beta, beta_null, ng_labels, 
                              weights=weights,
                              n_samples=n_samples_fit,
                              n_groups=n_groups,
                              degrees_of_freedom=degrees_of_freedom)
        
        colname_qval <- paste(contrast_name, "_qval", sep="")
        colname_pval <- paste(contrast_name, "_pval", sep="")
        
        results[colname_pval] <- pval
        qval <- stats::p.adjust(pval, method="BH")
        dim(qval) <- dim(pval)
        results[colname_qval] <- qval
    }
    
    return(results)
}

)
