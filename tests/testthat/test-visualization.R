library(moanin)
data(exampleData)
library(testthat)

context("moanin::visualization.R")

test_that("visualization::plot_splines_data", {
    moanin_model = create_moanin_model(data=testData,meta=testMeta)
    
    singleData = unlist(testData[1,])
    # Plot only one gene
    expect_silent(plot_splines_data(data=singleData , object=moanin_model))
    expect_silent(plot_splines_data(data=singleData , object=moanin_model, smooth=TRUE))
    expect_error(plot_splines_data(data=singleData , moanin_model,
				   smooth=TRUE, mfrow=c(1, 1),legend=TRUE))
    expect_silent(plot_splines_data(data=singleData ,  moanin_model,
				    smooth=TRUE, mfrow=c(1, 1),
				    legend=FALSE))

    # Plot several genes, using subset_data
    subset_data = row.names(testData)[2:5]
    expect_silent(plot_splines_data(data=testData, object=moanin_model,
				    subset_data=subset_data))

    # Don't provide data
    expect_silent(plot_splines_data(object=moanin_model,
                                    subset_data=subset_data))
    
    expect_error(plot_splines_data(data=testData, object=moanin_model,
				   subset_data=c(subset_data, "not_a_gene")))


    # Check that some options make sense
    expect_error(plot_splines_data(data=testData, moanin_model,
				   subset_data=subset_data,
				   mfrow=c(1, 1)))
		
    expect_error(plot_splines_data(data, moanin_model,
                                   subset_data=subset_data,
				   mfrow=c(1)))

})
