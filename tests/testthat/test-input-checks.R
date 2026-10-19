# Input that cannot work must stop a run at once, with a message that says what
# is wrong. These checks need no annotation, so they run before any long step.
# The tests of prepare_for_ORFquant() are in RiboseQC, which has its code now.

test_that("run_ORFquant() refuses wrong n_cores", {
  expect_error(run_ORFquant(for_ORFquant_file = "f", annotation_file = "a", n_cores = 0),
               "n_cores must be a number")
  expect_error(run_ORFquant(for_ORFquant_file = "f", annotation_file = "a"),
               "n_cores must be a number")
})

test_that("plot_orfquant_locus() warns when bam_files or col is given, because they are not used", {
  warned <- character()
  expect_error(withCallingHandlers(plot_orfquant_locus("TCN2", list(), "x.bam", col = "red"),
                                   warning = function(w) {
                                     warned <<- c(warned, conditionMessage(w))
                                     invokeRestart("muffleWarning")
                                   }),
               "old version|Please install")
  expect_equal(sub(":.*", "", warned), c("bam_files is not used", "col is not used"))
  expect_no_warning(expect_error(plot_orfquant_locus("TCN2", list()), "old version|Please install"))
})
