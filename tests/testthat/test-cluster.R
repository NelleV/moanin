data(exampleData)


context("moanin::cluster.R")

test_that("cluster::splines_kmeans", {
    random_seed = 42
    moanin_model = create_moanin_model(data=testData,meta=testMeta)
    expect_silent(moanin::splines_kmeans(moanin_model, n_init=1,
				         random_seed=random_seed))
})


test_that("cluster::splines_kmeans_score_and_label", {
    random_seed = 42
    moanin_model = create_moanin_model(data=testData,meta=testMeta)
    clustering_results = moanin::splines_kmeans(
	    moanin_model, n_init=1,
	    random_seed=random_seed, rescale=FALSE) #FIXME: temporary until figure out why creates NaN values

    expect_silent(scores_and_labels<-splines_kmeans_score_and_label(moanin_model, 
                               clustering_results))
    expect_equal(row.names(testData), row.names(scores_and_labels$scores))

    expect_silent(scores_and_labels2 <- splines_kmeans_score_and_label(moanin_model,
                   data=testData, clustering_results))
    expect_equal(scores_and_labels2, scores_and_labels)
    
    # Set a max score that we know is belove the max score found automatically
    max_score = max(scores_and_labels$scores) / 2
    scores_and_labels = splines_kmeans_score_and_label(
	    object=moanin_model,data=testData, clustering_results, max_score=(max_score))
    labels = scores_and_labels$labels
    scores = rowMin(scores_and_labels$scores[!is.na(labels), ])
    expect_true(max(scores) <= max_score)

    # Now, just do a ghost test with the rescale_separately_on

    expect_silent(
        splines_kmeans_score_and_label(
            moanin_model,
            clustering_results, rescale_separately=TRUE))

})

