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
test_that("moanin_model::discont_basis",{
    #test function
    x<-seq(0,10,length=100)
    basis<-discont_basis(x,discont_point=3, dfPre=3, dfPost=4, intercept=TRUE)
     
    # Use it in a moanin_model object instead of ns/bs:
    moanin <- create_moanin_model(data=testData, meta=testMeta,
        spline_formula=~Group:discont_basis(Timepoint,dfPre=3,
            dfPost=3,discont=20,intercept=TRUE)+0,
        df=6)
    
})