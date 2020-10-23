library(moanin)
library(limma)
data(exampleData)
context("moanin::de_analysis.R")

test_that("time-course DE analysis", {
    #change group variable to check that works everywhere
    colnames(testMeta)<-c("Index","condition","replicate","time")
    
    moanin_model = create_moanin_model(data=testData, 
        meta=testMeta,
        group_variable_name="condition",
        time_variable_name="time")
    # Create a second moanin model with counts
    moanin_model_counts = create_moanin_model(
        data=round(exp(testData)), meta=testMeta,
        log_transform=TRUE,
        group_variable_name="condition",
        time_variable_name="time")

    # Do only one contrast
    contrast_formulas = c("K-M")
    contrast = limma::makeContrasts(
	    contrasts=contrast_formulas,
	    levels=levels(droplevels(testMeta$condition)))

    expect_silent(moanin::DE_timecourse(
        moanin_model, contrasts=contrast,
	    use_voom_weights=FALSE))

    expect_silent(moanin::DE_timecourse(
        moanin_model_counts, contrasts=contrast,
	    use_voom_weights=TRUE))

    # Now test with several contrasts
    contrast_formulas = c("K-M", "C-M")
    contrasts = limma::makeContrasts(
	    contrasts=contrast_formulas,
	    levels=levels(droplevels(testMeta$condition)))

    expect_silent(moanin::DE_timecourse(moanin_model, contrasts=contrast,
	    use_voom_weights=FALSE))
        
    expect_silent(moanin::DE_timecourse(
        moanin_model_counts, contrasts=contrast,
	    use_voom_weights=TRUE))

    
    # Check that the names of the columns make sense.
    de_results = moanin::DE_timecourse(
	    moanin_model, contrasts=contrast,
	    use_voom_weights=FALSE)

    })
