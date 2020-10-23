library(moanin)
data(exampleData)
library(limma)

context("moanin::DE_timepoints.R")

test_that("DE_timepoints::create_timepoints_contrasts", {
    #change group variable to check that works everywhere
    colnames(testMeta)<-c("Index","condition","replicate","time")
    
    moanin_model = create_moanin_model(data=testData,meta=testMeta,
        group_variable_name="condition",
        time_variable_name="time")
    expect_silent(contrasts<-create_timepoints_contrasts( moanin_model,"M", "C"))

    # Now drop timepoint 6 in condition C. That should raise a warning
    mask = !((testMeta$condition == "C") & (testMeta$time == 6))
    msg = paste0("timepoint 6 is missing in condition C")
    expect_warning(create_timepoints_contrasts(moanin_model[,mask],
        "M", "C"), msg)
        
    expect_silent( DE_timepoints(moanin_model,
         contrasts=contrasts, use_voom_weights=FALSE))
    # Create a second moanin model with counts
    moanin_model_counts = create_moanin_model(
        data=round(exp(testData)), meta=testMeta,
        log_transform=TRUE,
        group_variable_name="condition",
        time_variable_name="time")
    expect_silent( DE_timepoints(moanin_model_counts,
          contrasts=contrasts, use_voom_weights=TRUE))
    
    #Check add replicate
    expect_silent( DE_timepoints(moanin_model_counts,
          contrasts=contrasts, use_voom_weights=TRUE,
          add_factors="replicate"))

    expect_error( DE_timepoints(moanin_model_counts,
        contrasts=contrasts, use_voom_weights=TRUE,
        add_factors="Reps"),"Error in creating the design matrix")

})
