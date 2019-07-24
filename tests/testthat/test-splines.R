library("moanin")

context("moanin::splines.R")

test_that("splines::fit_predict_splines", {
    data(shoemaker2015) 
    meta = shoemaker2015$meta
    data = shoemaker2015$data[1:5, ]

    splines_model = create_splines_model(meta)
    meta_prediction = create_meta_prediction(splines_model)

    expect_silent(fit_predict_splines(data, splines_model))
    expect_silent(fit_predict_splines(data, splines_model,
    				      meta_prediction=meta_prediction))

})

test_that("splines:align_data_onto_centroid", {
    set.seed(42)
    n_samples = 20
    n_genes = 5
    centroid = runif(n_samples)

    shift = runif(n_genes)
    scale = 1 + runif(n_genes)

    data = rep(centroid, each=n_genes)
    dim(data) = c(n_genes, n_samples)
    expect_equal(data, align_data_onto_centroid(data, centroid))

    # Ok, now let's make this a bit more complicated.
    scale = rep(scale, times=n_samples)
    dim(scale) = dim(data)
    shift = rep(shift, times=n_samples)
    dim(shift) = dim(data)

    shifted_scaled_data = scale * data + shift
    shifted_scaled_centroid = scale * centroid + shift
    expect_equal(data,
		 align_data_onto_centroid(shifted_scaled_data, centroid))
    expect_error(align_data_onto_centroid(data[, 1:10], centroid))
})

test_that("splines:score_genes_centroid", {
    set.seed(42)
    n_samples = 20
    n_genes = 5
    centroid = runif(n_samples)

    shift = runif(n_genes)
    scale = 1 + runif(n_genes)

    data = rep(centroid, each=n_genes)
    dim(data) = c(n_genes, n_samples)
    expect_silent(score_genes_centroid(data, centroid))
    expect_equal(sum(score_genes_centroid(data, centroid)), 0)

})

test_that("splines:rescale_values", {
    n_genes = 5
    data(shoemaker2015)
    data = shoemaker2015$data
    meta = shoemaker2015$meta
    data = data[1:n_genes,]

    rescaled_data = rescale_values(data, meta)
    expect_equal(rep(0, n_genes), as.vector(row_min(rescaled_data)))
    expect_equal(rep(1, n_genes), as.vector(row_max(rescaled_data)))
})
