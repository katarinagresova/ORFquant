# Input that cannot work must stop a run at once, with a message that says what
# is wrong. These checks need no annotation, so they run before any long step.

work <- tempfile("orfquant_input_")
dir.create(work)

test_that("prepare_for_ORFquant() and run_ORFquant() refuse wrong n_cores", {
  expect_error(prepare_for_ORFquant(annotation_file = "a", bam_file = "b", path_to_rl_cutoff_file = "c",
                                    n_cores = 0), "n_cores must be a number")
  expect_error(run_ORFquant(for_ORFquant_file = "f", annotation_file = "a", n_cores = 0),
               "n_cores must be a number")
  expect_error(run_ORFquant(for_ORFquant_file = "f", annotation_file = "a"),
               "n_cores must be a number")
})

test_that("prepare_for_ORFquant() names the files that do not exist", {
  msg <- tryCatch(prepare_for_ORFquant(annotation_file = file.path(work, "no_Rannot"),
                                       bam_file = file.path(work, "no.bam"),
                                       path_to_rl_cutoff_file = file.path(work, "no.tsv")),
                  error = function(e) conditionMessage(e))
  expect_match(msg, "don't exist")
  for (f in c("no_Rannot", "no.bam", "no.tsv")) expect_match(msg, f, fixed = TRUE)
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
