data(exampleData)
context("moanin::splines.R")

test_that("splines::fit_predict_splines", {

    moanin_model = create_moanin_model(data=testData,meta=testMeta)
    meta_prediction = create_meta_prediction(moanin_model)

    expect_silent(fit_predict_splines(testData, moanin_model))
    expect_silent(fit_predict_splines(testData, moanin_model,
    				      meta_prediction=meta_prediction))

})

test_that("splines:align_data_onto_centroid", {
    set.seed(42)
    n_samples = 20
    n_genes = 5
    centroid = runif(n_samples)

    ## don't do anything, just make data equal to centroid.
    ## (Shouldn't do anything)
    data = rep(centroid, each=n_genes)
    dim(data) = c(n_genes, n_samples)
    expect_equal(data, align_data_onto_centroid(data, centroid))
    expect_equal(data, align_data_onto_centroid(data, centroid, returnType="centroid"))
    
    # Ok, now let's make this a bit more complicated.
    shift_param = runif(n_genes)
    scale_param = 1 + runif(n_genes)
    scale_param = rep(scale_param, times=n_samples)
    dim(scale_param) = dim(data)
    shift_param = rep(shift_param, times=n_samples)
    dim(shift_param) = dim(data)

    shifted_scaled_data = scale_param * data + shift_param
    shifted_scaled_centroid = scale_param * centroid + shift_param
    expect_equal(data,
		 align_data_onto_centroid(shifted_scaled_data, centroid))
    expect_equal(shifted_scaled_data,
                 align_data_onto_centroid(shifted_scaled_data, centroid,returnType="centroid"))

    # Check error checking works
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
    moanin_model = moanin::create_moanin_model(data=testData,meta=testMeta)
    expect_silent(rescaled_data <- rescale_values(moanin_model,use_group=FALSE))
    expect_silent(rescaled_data <- rescale_values(moanin_model,data=moanin:::get_log_data(moanin_model),use_group=FALSE))
    expect_equal(rep(0, nrow(moanin_model)), as.vector(row_min(rescaled_data)))
    expect_equal(rep(1, nrow(moanin_model)), as.vector(row_max(rescaled_data)))

    #check different imputs
    expect_silent(rescaled_data2<- rescale_values(moanin_model,
                  data=get_log_data(moanin_model)[1:20,],use_group=FALSE))
    expect_equal(rescaled_data[1:20,],rescaled_data2)    
    expect_silent(rescaled_data3<- rescale_values(object=NULL,
                    data=assay(moanin_model)[1:20,]))
    expect_equal(rescaled_data3,rescaled_data2)    

})


test_that("splines::create_meta_prediction", {
    moanin_model = create_moanin_model(data=testData,meta=testMeta)
    expect_silent(create_meta_prediction(moanin_model))

    # Recreate moanin model without providing the formula
    # (will create warning in create_meta_prediction)
    basis = basis_matrix(moanin_model)
    moanin_model = create_moanin_model(data=testData,meta=testMeta, basis=basis)
    expect_warning(create_meta_prediction(moanin_model))
})
