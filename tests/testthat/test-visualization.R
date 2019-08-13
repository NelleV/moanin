library(moanin)
library(testthat)

context("moanin::visualization.R")

test_that("visualization::plot_centroid_individual", {
    data(shoemaker2015)
    data = as.vector(shoemaker2015$data[1, ])
    meta = shoemaker2015$meta

    splines_model = create_splines_model(meta)

    #expect_silent(plot_centroid_individual(data, splines_model))
    #expect_silent(plot_centroid_individual(data, splines_model, smooth=TRUE))

})


test_that("visualization::plot_splines_data", {
    data(shoemaker2015)
    data = as.matrix(shoemaker2015$data[5, ])
    meta = shoemaker2015$meta

    splines_model = create_splines_model(meta)

    expect_silent(plot_splines_data(data, splines_model))
    expect_silent(plot_splines_data(data, splines_model, smooth=TRUE))
    expect_error(plot_splines_data(data, splines_model,
				   smooth=TRUE, mfrow=c(1, 1)))
    expect_silent(plot_splines_data(data, splines_model,
				    smooth=TRUE, mfrow=c(1, 1),
				    legend=FALSE))

})
