library(moanin)
library(testthat)

context("moanin::visualization.R")

test_that("visualization::plot_centroid_individual", {
    data(shoemaker2015)
    data = as.vector(shoemaker2015$data[1, ])
    meta = shoemaker2015$meta

    moanin_model = create_moanin_model(meta)

    #expect_silent(plot_centroid_individual(data, moanin_model))
    #expect_silent(plot_centroid_individual(data, moanin_model, smooth=TRUE))

})


test_that("visualization::plot_splines_data", {
    data(shoemaker2015)
    data = as.matrix(shoemaker2015$data[5, ])
    meta = shoemaker2015$meta

    moanin_model = create_moanin_model(meta)

    # Plot only one gene
    expect_silent(plot_splines_data(data, moanin_model))
    expect_silent(plot_splines_data(data, moanin_model, smooth=TRUE))
    expect_error(plot_splines_data(data, moanin_model,
				   smooth=TRUE, mfrow=c(1, 1)))
    expect_silent(plot_splines_data(data, moanin_model,
				    smooth=TRUE, mfrow=c(1, 1),
				    legend=FALSE))

    # Plot several genes, using subset_data
    data = as.matrix(shoemaker2015$data[1:5,])
    subset_data = row.names(data)[2:5]
    expect_silent(plot_splines_data(data, moanin_model,
				    subset_data=subset_data))

    expect_error(plot_splines_data(data, moanin_model,
				   subset_data=c(subset_data, "not_a_gene")))


    # Check that some options make sense
    expect_error(plot_splines_data(data, moanin_model,
				   subset_data,
				   mfrow=c(1, 1)))
		
    expect_error(plot_splines_data(data, moanin_model,
				   subset_data,
				   mfrow=c(1)))

})