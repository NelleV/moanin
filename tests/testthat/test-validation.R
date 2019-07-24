library(moanin)

context("moanin::validation")

test_that("validation:check_meta", {
    data(shoemaker2015)
    meta = shoemaker2015$meta
    expect_silent(check_meta(meta))

    meta_without_group = subset(meta, select=c("Replicate", "Timepoint"))
    meta_without_time = subset(meta, select=c("Group", "Replicate"))
    expect_error(check_meta(meta_without_group))
    expect_error(check_meta(meta_without_time))

    meta_without_replicate = subset(meta, select=c("Timepoint", "Group"))
    expect_silent(check_meta(meta_without_replicate))
    expect_error(check_meta(meta_without_replicate, check_replicates=TRUE))

    meta_timepoint_not_numeric = meta
    meta_timepoint_not_numeric$Timepoint = as.factor(meta_timepoint_not_numeric$Timepoint)
    expect_error(check_meta(meta_timepoint_not_numeric))

    meta_group_not_factor = meta
    meta_group_not_factor$Group = as.character(meta$Group)
    expect_error(check_meta(meta_group_not_factor))


})

test_that("validation:check_data_meta", {
    data(shoemaker2015)
    data = shoemaker2015$data
    meta = shoemaker2015$meta

    data = data[1:10, ]

    # Running check_data_meta should work fine
    expect_silent(check_data_meta(data, meta))

    expect_error(
	check_data_meta(data[, 2:10], meta))
})


test_that("validation:check_is_2d", {
    data_2d = matrix(1:9, nrow=3)
    expect_silent(check_is_2d(data_2d))
    expect_silent(check_is_2d(as.data.frame(data_2d)))
    expect_error(check_is_2d(1:9))
})
