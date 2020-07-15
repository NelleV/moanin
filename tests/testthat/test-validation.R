library(moanin)
data(exampleData)
context("moanin::validation")


test_that("validation:check_data_meta", {
    moanin_model = create_moanin_model(data=testData,meta=testMeta)
    
    # Running check_data_meta should work fine
    expect_silent(check_data_meta(testData, moanin_model))

    expect_error(
	check_data_meta(testData[, 2:10], moanin_model))
})


test_that("validation:check_is_2d", {
    data_2d = matrix(1:9, nrow=3)
    expect_silent(check_is_2d(data_2d))
    expect_silent(check_is_2d(as.data.frame(data_2d)))
    expect_error(check_is_2d(1:9))
})


test_that("validation:is_contrasts", {
    contrasts_formula = c("M-K", "M-C")
    moanin_model = create_moanin_model(data=testData,meta=testMeta)
    expect_silent(contrasts <-
	    is_contrasts(contrasts_formula, moanin_model))
    expect_silent(is_contrasts(contrasts,  moanin_model))

    contrasts = limma::makeContrasts(
	    contrasts=contrasts_formula, levels=levels(testMeta$Group))
    expect_silent(
	    is_contrasts(contrasts,  moanin_model))

})
