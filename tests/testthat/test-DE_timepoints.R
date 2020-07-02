library(moanin)
library(timecoursedata)
library(limma)

context("moanin::DE_timepoints.R")

test_that("DE_timepoints::create_timepoints_contrasts", {
    data(shoemaker2015)
    data = shoemaker2015$data
    meta = shoemaker2015$meta
    
    moanin_model = create_moanin_model(meta)
    expect_silent(create_timepoints_contrasts("M", "C", moanin_model))

    # Now drop timepoint 6 in condition C. That should raise a warning
    mask = !((meta$Group == "C") & (meta$Timepoint == 6))
    moanin_model = create_moanin_model(meta[mask,])
    msg = paste0("moanin::create_timepoints_contrasts: timepoint",
		 " 6 is missing in condition C")
    expect_warning(create_timepoints_contrasts("M", "C", moanin_model), msg)

})
