test_that("check_requirements works", {

  run_dir <- paste0(tempdir(), do.call(paste0, replicate(4, sample(LETTERS, 1, TRUE), FALSE)),"chk")
  if (!dir.exists(run_dir)) {
    dir.create(run_dir)
  }
  models_dir <- paste0(tempdir(), do.call(paste0, replicate(4, sample(LETTERS, 1, TRUE), FALSE)),"chk")
  if (!dir.exists(models_dir)) {
    dir.create(models_dir)
  }
  # copy nmfe74.bat to run_dir

  # will get model files from copy files
  source_nmfe <- file.path(system.file(package = "mbbe", "test_files", "check_requirements","nmfetest.bat"))

  source_models <- system.file(package = "mbbe", "test_files", "copy_files")
  #$DATA  U:\\fda\\mbbe\\inst\\examples\\data_seq.csv does not exist so we must
  # update path in control file and write .mod files with new data path
  model1 <- suppressWarnings(readLines(file.path(source_models,"model1.mod")))
  model2 <- suppressWarnings(readLines(file.path(source_models,"model2.mod")))
  data_file <- system.file(package = "mbbe", "examples", "data_seq.csv")
  file.copy(from = data_file, to = models_dir, overwrite = TRUE)
  data_file <- gsub("\\\\", "/", file.path(models_dir, "data_seq.csv"))
  model1[5] <- sub("\\$DATA[[:space:]]*(.*?)[[:space:]]*IGNORE=@",
                   paste0("$DATA ", data_file, " IGNORE=@"),
                   model1[5])
  model2[5] <- sub("\\$DATA[[:space:]]*(.*?)[[:space:]]*IGNORE=@",
                   paste0("$DATA ", data_file, " IGNORE=@"),
                   model2[5])
  model1_temp_path <- file.path(models_dir,"model1.mod")
  model2_temp_path <- file.path(models_dir,"model2.mod")
  writeLines(model1, model1_temp_path)
  writeLines(model2, model2_temp_path)

  model_list <-  c(model1_temp_path, model2_temp_path)
  dir.create(file.path(run_dir,"MBBEsim1"))
  source_sim_data <- file.path(system.file(package = "mbbe", "test_files", "get_parms","data_sim.csv"))


  with_mock(
    askYesNo = function(...) {
      return(TRUE)
    },
    test_ReturnedValue <-  check_requirements(run_dir, 1,
                                              model_list,
                                              4,
                                              c(1,2),
                                              c(3,4),
                                              source_nmfe,
                                              TRUE,
                                              source_sim_data,
                                              FALSE)
  )
  expect_true(test_ReturnedValue$rval)
  unlink(run_dir, recursive = TRUE)
  unlink(models_dir,  recursive = TRUE)
})
