library("moanin")

context("moanin::cluster_validation.R")

test_that("cluster_validation::consensus_matrix", {
    set.seed(42)
    n_genes = 40
    n_clusters = 5
    n_labels = 20
    labels = sample.int(n_labels, n_clusters * n_genes, replace=TRUE)
    dim(labels) = c(n_genes, n_clusters)
    colnames(labels) = lapply(1:n_clusters, function(x){paste0("C", x)})
    row.names(labels) = lapply(1:n_genes, function(x){paste0("G", x)})

    expect_silent(consensus_matrix(labels))
})
