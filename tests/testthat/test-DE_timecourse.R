library(moanin)
library(limma)
data(exampleData)
context("moanin::de_analysis.R")

test_that("time-course DE analysis", {

    moanin_model = create_moanin_model(data=testData,meta=testMeta)

    # Do only one contrast
    contrast_formulas = c("K-M")
    contrast = limma::makeContrasts(
	contrasts=contrast_formulas,
	    levels=levels(droplevels(testMeta$Group)))

    expect_silent(moanin::DE_timecourse(moanin_model, contrasts=contrast,
	    use_voom_weights=FALSE))

    # Now test with several contrasts
    contrast_formulas = c("K-M", "C-M")
    contrasts = limma::makeContrasts(
	    contrasts=contrast_formulas,
	    levels=levels(droplevels(testMeta$Group)))

    expect_silent(moanin::DE_timecourse(moanin_model, contrasts=contrast,
	    use_voom_weights=FALSE))
    
    # Check that the names of the columns make sense.
    de_results = moanin::DE_timecourse(
	    moanin_model, contrasts=contrast,
	    use_voom_weights=FALSE)

})
