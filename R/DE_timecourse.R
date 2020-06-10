# This file contains EDGE code adapted to work with limma contrasts

center_data = function(y, ng_labels){
  for(g in levels(ng_labels)){
    whKeep = which(ng_labels == g)
    sub_mean = rowMeans(y[, whKeep])
    y[, whKeep] = y[, whKeep] - sub_mean
  }
  return(y)
}

compute_beta_null = function(basis, beta, contrasts_coef){
  ng = length(contrasts_coef)
  df = ncol(basis) / ng
  contrasts_coef_ = rep(contrasts_coef, times=df)

  # Reshape b so that each row corresponds to a group and drop the intercept
  b_ = array(beta, dim=c(dim(beta)[1], ng, df))
  # First start by constructing the matrix T
  t_ = apply(b_, 1, function(x){colSums(x*contrasts_coef)})

  # The observations can not be assumed to be balanced...
  # We need to get rid of the intercept for this part
  # FIbasisME don't invert this matrix…
  part_K = MASS::ginv(t(basis) %*% basis)
  K = part_K * contrasts_coef_**2

  # We now need to sum all elements associated to the same pairs of splines.
  # Which is, in a particular case every four elements in both directions.

  K = sapply(1:df, function(jg) rowSums(K[, (jg-1)*ng + 1:ng]))
  K = sapply(1:df, function(jg) colSums(K[(jg-1)*ng + 1:ng,]))

  T_ = MASS::ginv(K) %*% t_

  # We got T. Now, let's move on to the rest
  tmp = as.array(rep(as.vector(T_), each=ng), dim=c(1, 1, 1))
  dim(tmp) = c(ng * df, dim(beta)[1])
  C_ = contrasts_coef * part_K %*% tmp

  C_ = t(C_)
  dim(C_) = c(dim(beta)[1], ng, df)
  beta_null = b_ - C_

  # Last step: reshape beta null so that it is of the same shape as beta
  dim(beta_null) = dim(beta)
  return(beta_null)
}


lrtStat = function(resNull, resFull, ng_labels=NULL) {
  # FIbasisME I'm pretty sure that in the case of contrasts, the degrees of
  # freedom computed here are wrong as they include part of the data that is
  # not used for the test. This needs to be fixed
  stat = 0
  if(is.null(ng_labels)){
    ss0 = rowSums(resNull^2)
    ss1 = rowSums(resFull^2)
    n = ncol(resNull)
    stat = stat + n * (ss0 - ss1)/(ss1)

  }else{
    for(g in levels(ng_labels)){
	whKeep = which(ng_labels == g)
	sub_resNull = resNull[,whKeep]
	sub_resFull = resFull[,whKeep]

	# Somehow the two lines above don't return the same object depending on
	# the dimension of resNull and resFull, so need to distinguish the case
	# where there is only one observation in data.
	if(is.null(dim(sub_resNull))){
	    ss0 = sum(sub_resNull^2)
	    ss1 = sum(sub_resFull^2)
	    n = length(sub_resNull)
	}else{
	    ss0 = rowSums(sub_resNull^2)
	    ss1 = rowSums(sub_resFull^2)
	    n = ncol(sub_resNull)
	}
	stat = stat + n * (ss0 - ss1)/(ss1)
    }
  }

  return(stat)
}

compute_pvalue = function(basis, y, beta, beta_null, ng_labels,
			  n_groups=NULL,
			  n_samples=NULL,
			  degrees_of_freedom=NULL,
			  statistics="lrt",
			  df2=NULL, weights=NULL){

    fitFull = beta %*% t(basis)

    fitNull = beta_null %*% t(basis)

    if(!is.null(weights)){
      resNull = weights^(1/2) * (y - fitNull)
      resFull = weights^(1/2) * (y - fitFull)
    }else{
      resNull = y - fitNull
      resFull = y - fitFull
    }

    # estimate degrees of freedom.
    if(is.null(n_groups)){
        n_groups = nlevels(ng_labels)
	# FIbasisME Raise warning
    }

    if(is.null(n_samples)){
	n_samples = ncol(basis)
	# FIbasisME raise warning
    }
    
    df = degrees_of_freedom

    if(statistics == "ftest"){
	stat = lrtStat(resNull, resFull)
	if(is.null(df2)){
	    # FIbasisME Check this.
	    df2 = n_samples - degrees_of_freedom * n_groups
	}
	df1 = df
	pval = stats::pf(stat * df2 / df1, df1=df1, df2=df2, lower.tail=FALSE)
    }else{
	lstat = lrtStat(resNull, resFull, ng_labels=ng_labels)
	pval = stats::pchisq(lstat, df=degrees_of_freedom, lower.tail=FALSE)
    }
    return(pval)
}

summarise = function(basis, ng_levels) {
  basis_mean = matrix(,nrow=nrow(basis), ncol=nlevels(ng_levels))
  colnames(basis_mean) = levels(ng_levels)
  rownames(basis_mean) = rownames(basis)

  for(g in levels(ng_levels)){
    whKeep = which(ng_levels == g)
    if(length(whKeep) > 1){
      basis_mean[, g] = rowMeans(basis[, whKeep])
    }else if(length(whKeep) != 0){
      basis_mean[, g] = basis[, whKeep]
    }
  }
  return(basis_mean)
}

#' Run spline models with contrasts.
#' 
#' @param data The data matrix, where genes (features) are on the rows and
#'   samples on the columns.
#' @param moanin_model object of class \code{moanin_model}, an object containing
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
#'	Whether to use voom weights.
#' @seealso \code{\link[limma]{makeContrasts}}, \code{\link{create_moanin_model}}, \code{\link{DE_timepoints}}, \code{\link[edge]{edge}}
#' @return A \code{data.frame} with two columns for each of the contrasts given
#'   in \code{contrasts}, corresponding to the raw p-value of the contrast for
#'   that gene (\code{_pval}) and the adjusted p-value (\code{_qval}). The
#'   adjusted p-values are FDR-adjusted based on the Benjamini-Hochberg method,
#'   as implemented in \code{\link[stats]{p.adjust}}. The adjustment is done
#'   across all p-values for all contrasts calculated.
#' @details The implementation of the spline fit and the calculation of p-values
#'   was based on code from \code{\link[edge]{edge}}, and expanded to
#'   enable handling of comparisons of groups via contrasts.
#' @export
DE_timecourse = function(data, moanin_model,
			 contrasts,
			 center=FALSE,
			 use_voom_weights=TRUE){

    basis = moanin_model$basis
    meta = moanin_model$meta

    ng = nlevels(meta$Group)
    ng_labels = meta$Group

    check_data_meta(data, meta)
    data = as.matrix(data)

    contrasts = is_contrasts(contrasts, meta)

    if(use_voom_weights){
        y = edgeR::DGEList(counts=data)
        y = edgeR::calcNormFactors(y, method="upperquartile")
        v = limma::voom(y, contrasts, plot=FALSE)
        weights = limma::lmFit(v)
    }else{
        weights = NULL
    }

    y = data

    if(center){
        y = center_data(y)
        basis = t(center_data(t(basis)))
    }
    
    beta = fit_splines(y, moanin_model, weights=weights)
    
    if(dim(contrasts)[1] != ng){
        stop("The contrast coef vector should be of the same size" +
                 " as the number of groups")
    }
    
    results = data.frame(row.names=row.names(data))
    for(col in 1:ncol(contrasts)){
        contrast = contrasts[, col]
        
        # Create the name of the column
        contrast_name = colnames(contrasts)[col]
        contrast_name = gsub(" ", "", contrast_name, fixed=TRUE)
        
        # Get the number of samples used for this particular contrast:
        groups_of_interest = names(contrast)[contrast != 0]
        n_samples_fit = sum(with(meta, Group %in% groups_of_interest))
        n_groups = length(groups_of_interest)
        degrees_of_freedom = dim(basis)[2] / ng
        
        beta_null = compute_beta_null(basis, beta, contrast)
        
        pval = compute_pvalue(basis, y, beta, beta_null, ng_labels, weights=weights,
                              n_samples=n_samples_fit,
                              n_groups=n_groups,
                              degrees_of_freedom=degrees_of_freedom)
        
        colname_qval = paste(contrast_name, "_qval", sep="")
        colname_pval = paste(contrast_name, "_pval", sep="")
        
        results[colname_pval] = pval
        qval = stats::p.adjust(pval, method="BH")
        dim(qval) = dim(pval)
        results[colname_qval] = qval
    }
    
    return(results)
}


