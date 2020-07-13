library(moanin)
data(exampleData)
library(limma)

context("moanin::DE_timepoints.R")

test_that("DE_timepoints::create_timepoints_contrasts", {

    moanin_model = create_moanin_model(data=testData,meta=testMeta)
    expect_silent(create_timepoints_contrasts( moanin_model,"M", "C"))

    # Now drop timepoint 6 in condition C. That should raise a warning
    mask = !((testMeta$Group == "C") & (testMeta$Timepoint == 6))
    msg = paste0("timepoint 6 is missing in condition C")
    expect_warning(create_timepoints_contrasts(moanin_model[,mask],"M", "C"), msg)

})
