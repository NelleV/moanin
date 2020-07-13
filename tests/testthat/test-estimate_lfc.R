library("moanin")

context("moanin::de_analysis_utils.R")
data(exampleData)

test_that("Estimating log fold change smoke tests", {
    # This is just a smoke test.
    moanin_model = moanin::create_moanin_model(data=testData,meta=testMeta)

    # Reduce the data set
    methods = eval(formals(estimate_log_fold_change)$method)

    contrast_formula = c("C-K")
    contrasts = limma::makeContrasts(contrasts=contrast_formula,
				     levels=droplevels(testMeta$Group))
    for(method in methods){
        expect_silent(
	        estimate_log_fold_change(moanin_model, contrasts, method=method))
	    expect_silent(
	        estimate_log_fold_change(moanin_model, contrast_formula, method=method))
    }

    # Now same test, but several contrasts
    contrasts = limma::makeContrasts(contrasts=c("C-K", "C-M"),
				     levels=droplevels(testMeta$Group))
    for(method in methods){
        expect_silent(
	    estimate_log_fold_change(
		    object=moanin_model, contrasts, method=method))
    }

})

test_that("Estimating log fold change with unknown error", {

    moanin_model = moanin::create_moanin_model(data=testData,meta=testMeta)
    
    # Reduce the data set
    contrasts = c("C-K")
    expect_error(estimate_log_fold_change(moanin_model,
					  contrast, method="hahaha"))
})


test_that("Estimating log fold change", {
    data<-testData
    data[, testMeta$Group == "C"] = 0
    data[, testMeta$Group == "K"] = 1
    moanin_model = moanin::create_moanin_model(data=data,meta=testMeta)
    # Reduce the data set
    contrasts = limma::makeContrasts(contrasts=c("C-K"), levels=testMeta$Group)
    lfc_max = estimate_log_fold_change(moanin_model, contrasts=contrasts, method="max")
    expect_equal(1, max(lfc_max))
    expect_equal(1, min(lfc_max))

})
