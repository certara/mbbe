test_that("copy_model_files works", {

  run_dir <- paste0(tempdir(), do.call(paste0, replicate(5, sample(LETTERS, 1, TRUE), FALSE)))
  if (!dir.exists(run_dir)) {
    dir.create(run_dir)
  }
  #when package is installed files in inst go to r library/package

  model1_file <- test_path("../test_files", "copy_files","model1.mod")
  model2_file <- test_path("../test_files", "copy_files","model2.mod")
  # make list of model files
  model_files <- c(model1_file, model2_file)
  nmodels <- mbbe:::copy_model_files(model_files, run_dir)

  con <- file(model1_file, "r")
  suppressWarnings(source_file1 <- readLines(con, encoding = "UTF-8"))
  close(con)
  con <- file(model2_file, "r")
  suppressWarnings(source_file2 <- readLines(con, encoding = "UTF-8"))
  close(con)
  con <- file(file.path(run_dir, "model1", "bs1.mod"), "r")
  suppressWarnings(run_file1 <- readLines(con, encoding = "UTF-8"))
  close(con)
  con <- file(file.path(run_dir,"model2", "bs2.mod"), "r")
  suppressWarnings(run_file2 <- readLines(con, encoding = "UTF-8"))
  close(con)
  expect_equal(nmodels,2)
  expect_equal(source_file1, run_file1)
  expect_equal(source_file2, run_file2)

  unlink(run_dir, recursive = TRUE)
})
