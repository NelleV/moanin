library(moanin)

context("moanin::utils.R")


test_that("Testing rowMin & rowMax", {
    lfc = matrix(c(0, 0, 2, -1, 0, 0), nrow=2, ncol=3)
    expect_equal(c(2, 0),
		 rowMax(lfc))
    expect_equal(c(0, -1), rowMin(lfc))
    expect_error(rowMin(c(0, -1)))
    expect_error(rowMax(c(0, -1)))
})
