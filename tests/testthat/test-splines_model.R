library("moanin")
data("exampleData")

context("moanin::moanin_model.R")

test_that("moanin_model::create_moanin_model", {
    expect_silent(create_moanin_model(data=testData,meta=testMeta))
    formula = ~Group:splines::ns(Timepoint) + 0
    expect_silent(create_moanin_model(data=testData,meta=testMeta, 
                                      spline_formula=formula))
    basis = stats::model.matrix(formula, data=testMeta)
    expect_silent(create_moanin_model(data=testData,meta=testMeta, 
                                      basis_matrix=basis))

    expect_error(create_moanin_model(data=testData,meta=testMeta, 
                       basis_matrix=basis, spline_formula=formula))
})
