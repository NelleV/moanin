library("moanin")

context("moanin::de_analysis_utils.R")


test_that("Estimating log fold change smoke tests", {
    # This is just a smoke test.
    data(shoemaker2015)
    data = shoemaker2015$data
    meta = shoemaker2015$meta

    splines_model = moanin::create_splines_model(meta)

    # Reduce the data set
    data = data[1:10, ]
    methods = moanin:::ALL_LFC_METHODS

    contrasts = limma::makeContrasts(contrasts="C-K",
				     levels=meta$Group)
    for(method in methods){
        expect_silent(
	    estimate_log_fold_change(
		data, splines_model, contrasts, method=method))
    }

    # Now same test, but several contrasts
    contrasts = limma::makeContrasts(contrasts=c("C-K", "C-M"),
				     levels=meta$Group)
    for(method in methods){
        expect_silent(
	    estimate_log_fold_change(
		data, splines_model, contrasts, method=method))
    }

})

test_that("Estimating log fold change with unknown error", {
    data(shoemaker2015)
    data = shoemaker2015$data
    meta = shoemaker2015$meta

    splines_model = moanin::create_splines_model(meta)

    # Reduce the data set
    data = data[1:10, ]
    contrasts = c("C-K")
    expect_error(estimate_log_fold_change(data, splines_model,
					  contrast, method="hahaha"))
})


test_that("Estimating log fold change", {
    data(shoemaker2015)
    data = shoemaker2015$data
    meta = shoemaker2015$meta

    splines_model = moanin::create_splines_model(meta)

    # Reduce the data set
    data = data[1:10, ]
    contrasts = limma::makeContrasts(contrasts=c("C-K"), levels=meta$Group)
    data[, meta$Group == "C"] = 0
    data[, meta$Group == "K"] = 1
    lfc_max = estimate_log_fold_change(data, splines_model, contrasts=contrasts, method="max")
    expect_equal(1, max(lfc_max))
    expect_equal(1, min(lfc_max))

})
