library(moanin)
library(limma)
library(timecoursedata)
data(shoemaker2015)

context("moanin::de_analysis.R")

test_that("time-course DE analysis", {
    data = shoemaker2015$data
    meta = shoemaker2015$meta

    moanin_model = create_moanin_model(meta)

    # Limit ourselves to the top 5000 genes
    data = data[1:5000, ]

    # Do only one contrast
    contrast_formulas = c("K-M")
    contrast = limma::makeContrasts(
	contrasts=contrast_formulas,
	levels=levels(meta$Group))

    expect_silent(moanin::DE_timecourse(
	as.matrix(data), moanin_model, contrasts=contrast,
	use_voom_weights=FALSE))

    # Now test with several contrasts
    contrast_formulas = c("K-M", "C-M")
    contrasts = limma::makeContrasts(
	contrasts=contrast_formulas,
	levels=levels(meta$Group))

    expect_silent(moanin::DE_timecourse(
    	as.matrix(data), moanin_model, contrasts=contrast,
	use_voom_weights=FALSE))
    
    # Check that the names of the columns make sense.
    de_results = moanin::DE_timecourse(
	as.matrix(data), moanin_model, contrasts=contrast,
	use_voom_weights=FALSE)

})
