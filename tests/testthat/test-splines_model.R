library("moanin")

context("moanin::moanin_model.R")

test_that("moanin_model::create_moanin_model", {
    data(shoemaker2015)
    meta = shoemaker2015$meta
    expect_silent(create_moanin_model(meta))
    formula = ~Group:splines::ns(Timepoint) + 0
    expect_silent(create_moanin_model(meta, formula=formula))
    basis = stats::model.matrix(formula, data=meta)
    expect_silent(create_moanin_model(meta, basis=basis))

    expect_error(create_moanin_model(meta, basis=basis, formula=formula))
})
