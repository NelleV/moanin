library(moanin)
library(testthat)


context("moanin::cluster.R")

test_that("cluster::splines_kmeans", {
    data(shoemaker2015)
    data = shoemaker2015$data
    meta = shoemaker2015$meta
    random_seed = 42

    data = data[1:500, ]
    moanin_model = create_moanin_model(meta)
    expect_silent(moanin::splines_kmeans(data, moanin_model, n_init=1,
				         random_seed=random_seed))
})


test_that("cluster::splines_kmeans_score_and_label", {
    data(shoemaker2015)
    data = shoemaker2015$data
    meta = shoemaker2015$meta
    random_seed = 42

    data = data[1:500, ]
    moanin_model = create_moanin_model(meta)
    clustering_results = moanin::splines_kmeans(
	    data, moanin_model, n_init=1,
	    random_seed=random_seed)

    expect_silent(splines_kmeans_score_and_label(data, clustering_results))
    scores_and_labels = splines_kmeans_score_and_label(data, clustering_results)
    expect_equal(row.names(data), row.names(scores_and_labels$scores))

    # Set a max score that we know is belove the max score found automatically
    max_score = max(scores_and_labels$scores) / 2
    scores_and_labels = splines_kmeans_score_and_label(
	data, clustering_results, max_score=(max_score))
    labels = scores_and_labels$labels
    scores = rowMin(scores_and_labels$scores[!is.na(labels), ])
    expect_true(max(scores) <= max_score)

    # Now, just do a ghost test with the rescale_separately_on

    expect_silent(
        splines_kmeans_score_and_label(
            data,
            clustering_results, rescale_separately_on="Group"))

})

