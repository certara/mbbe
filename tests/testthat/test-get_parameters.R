test_that("get_parameters works", {

  run_dir <- paste0(tempdir(), do.call(paste0, replicate(5, sample(LETTERS, 1, TRUE), FALSE)))
  if (!dir.exists(run_dir)) {
    dir.create(run_dir)
  }
  # copy everything to run_dir

  souce_dir <- system.file(package = "mbbe", "test_files", "get_parms")
  dir.create(file.path(run_dir,"model1"))
  dir.create(file.path(run_dir,"model1","1"))
  file.copy(file.path(system.file(package = "mbbe", "test_files", "get_parms", "model1", "1"),"bsSamp1_1.lst"),
            file.path(run_dir,"model1","1"))
  file.copy(file.path(system.file(package = "mbbe", "test_files", "get_parms", "model1", "1"),"bsSamp1_1.ext"),
            file.path(run_dir,"model1","1"))
  file.copy(file.path(system.file(package = "mbbe", "test_files", "get_parms", "model1", "1"),"bsSamp1_1.xml"),
            file.path(run_dir,"model1","1"))
  dir.create(file.path(run_dir,"model2"))
  dir.create(file.path(run_dir,"model2","1"))
  file.copy(file.path(system.file(package = "mbbe", "test_files", "get_parms", "model2", "1"),"bsSamp2_1.lst"),
            file.path(run_dir,"model2","1"))
  file.copy(file.path(system.file(package = "mbbe", "test_files", "get_parms", "model2", "1"),"bsSamp2_1.ext"),
            file.path(run_dir,"model2","1"))
  file.copy(file.path(system.file(package = "mbbe", "test_files", "get_parms", "model2", "1"),"bsSamp2_1.xml"),
            file.path(run_dir,"model2","1"))


  referenceBICs <- read.csv(file.path(souce_dir, "BICS.csv"))
  get_parameters(run_dir, 2, 1, 0.2, 999999, TRUE)

  testBICs <- read.csv(file.path(run_dir,"BICS.csv"))
  expect_equal(testBICs, referenceBICs)

})
