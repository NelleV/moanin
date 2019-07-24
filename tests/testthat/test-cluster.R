library(moanin)
library(testthat)


context("moanin::cluster.R")

test_that("cluster::splines_kmeans", {
    data(shoemaker2015)
    data = shoemaker2015$data
    meta = shoemaker2015$meta
    random_seed = 42

    data = data[1:500, ]
    splines_model = create_splines_model(meta)
    expect_silent(moanin::splines_kmeans(data, splines_model, n_init=1,
				         random_seed=random_seed))
})


test_that("cluster::splines_kmeans_score_and_label", {
    data(shoemaker2015)
    data = shoemaker2015$data
    meta = shoemaker2015$meta
    random_seed = 42

    data = data[1:500, ]
    splines_model = create_splines_model(meta)
    clustering_results = moanin::splines_kmeans(
	data, splines_model, n_init=1,
	random_seed=random_seed)

    expect_silent(splines_kmeans_score_and_label(data, clustering_results))
    scores_and_labels = splines_kmeans_score_and_label(data, clustering_results)
    expect_equal(row.names(data), row.names(scores_and_labels$scores))
})

