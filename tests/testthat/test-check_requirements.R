test_that("check_requirements works", {

  run_dir <- paste0(tempdir(), do.call(paste0, replicate(4, sample(LETTERS, 1, TRUE), FALSE)),"chk")
  if (!dir.exists(run_dir)) {
    dir.create(run_dir)
  }
  # copy nmfe74.bat to run_dir

  # will get model files from copy files
  source_nmfe <- file.path(system.file(package = "mbbe", "test_files", "check_requirements","nmfetest.bat"))

  source_models <- system.file(package = "mbbe", "test_files", "copy_files")

  model1 <- file.path(source_models,"model1.mod")
  model2 <- file.path(source_models,"model2.mod")
  model_list <-  c(model1, model2)
  dir.create(file.path(run_dir,"MBBEsim1"))
  source_sim_data <- file.path(system.file(package = "mbbe", "test_files", "get_parms","data_sim.csv"))
  Reference_ReturnedValue <- list()
  Reference_ReturnedValue$rval <- TRUE
  Reference_ReturnedValue$msg <- "Passes requirements check"

  test_ReturnedValue <-  check_requirements(run_dir, 1,
                                            model_list,
                                            4,
                                            c(1,2),
                                            c(3,4),
                                            source_nmfe,
                                            TRUE,
                                            source_sim_data,
                                            FALSE)

  expect_equal(Reference_ReturnedValue, test_ReturnedValue)

})
