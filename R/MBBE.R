#' @importFrom stats coef confint filter lag lm na.exclude qchisq reshape runif qf
#' @importFrom utils read.csv read.table write.csv capture.output head tail
#' @importFrom magrittr %>%
#' @rawNamespace import(future, except = run)
#' @importFrom processx run
#' @import stringr
#' @import ps
NULL

`%>%` <- magrittr::`%>%`

#' delete_files
#'
#' @param folder A folder from which to delete non-essential NONMEM files
#' @details Files to be deleted are:
#' FDATA",
#' FCON,
#' FSUBS,
#' FSUBS.o,
#' FSIZES,
#' FREPORT,
#' FSTREAM,
#' GFCOMPILE.BAT,
#' INTER,
#' nmprd4p.mod,
#' and file with the following extensions:
#' exe, f90, grd, shk, cpu, shm, lnk, phi
#'
#'@export
delete_files <- function(folder) {
  # need .ext for use_identifiable, keep .lst, .xml, .mod, FMSG
  # but can't delete all, need to keep $TABLE
  try({
    delfiles <-
      c(
        "FDATA",
        "FCON",
        "FSUBS",
        "FSUBS.o",
        "FSIZES",
        "FREPORT",
        "FSTREAM",
        "GFCOMPILE.BAT",
        "INTER",
        "nmprd4p.mod"
      )

    ExtensionsToDelete <- c("exe", "f90", "grd", "shk", "cpu", "shm", "lnk", "phi")
    ExtensionsToDelete <- paste0("(*.", ExtensionsToDelete, ")")
    ExtensionsToDeleteCollapsed <- paste(ExtensionsToDelete, collapse = "|")
    delfiles <-
      c(
        delfiles,
        dir(path = folder, pattern = ExtensionsToDeleteCollapsed)
      )

    file.remove(file.path(folder, delfiles))
    if (dir.exists(file.path(folder, "temp_dir"))) {
      unlink(file.path(folder, "temp_dir"),
        recursive = TRUE,
        force = TRUE
      )
    }
  })
}

#' check_requirements
#'
#' Check if all the requirements for model based BE assessment are present
#'
#' @param run.dir source  model files
#' @param samp_size number of samples - must be integer
#' @param model_list list of model control files to be used for model averaging
#' @param ngroups number of groups for BE testing in data set
#' @param reference_groups Which groups are reference
#' @param test_groups Which groups are Test
#' @param nmfe_path nmfe??.bat
#' @param use_check_identifiable - logical, whether to check for identifability, will check if SADDLE_RESET is in $EST
#' @param use_simulation_data logical, if the simulation will be done with a different data set than the bootstrap
#' @param simulation_data_path if use_simulation_data, this is the path to the data file
#'
#' @details Check to see if:
#' 1. Does that .mod file contain ;;;;.*Start EST
#' 2. Is the sum of reference_groups and test_groups == ngroups?
#' 3. Any duplicates in Reference and Test groups?
#' 4. check if data file is present
#' 5. check if nmfe path is correct
#' 6. check for saddle_reset if requested for check_identifiability
#' 7. check for repeat IDs in data set
#' 8. if use_simulation_data, see if data available
check_requirements <- function(run_dir,
                               samp_size,
                               model_list,
                               ngroups,
                               reference_groups,
                               test_groups,
                               nmfe_path,
                               use_check_identifiable,
                               use_simulation_data,
                               simulation_data_path = NULL) {
  msg <- remove_old_files(run_dir, samp_size)

  if (!file.exists(nmfe_path)) {
    msg <-
      paste(msg,
        paste0("Cannot find nmfe?? at ", nmfe_path),
        sep = "\n"
      )
  }

  # check number in Reference and Test groups
  if (sum(length(reference_groups), length(test_groups)) != ngroups) {
    msg <-
      paste(
        msg,
        paste(
          "Number of Reference groups",
          length(reference_groups),
          "+ Test groups",
          length(test_groups),
          "doesn't equal the number of groups",
          ngroups
        ),
        sep = "\n"
      )
  }

  # no duplicated in Reference and Test groups
  if (anyDuplicated(c(reference_groups, test_groups)) > 0) {
    msg <-
      paste(msg,
        "There are duplicated group numbers between Reference and Test group",
        sep = "\n"
      )
  }

  if (anyDuplicated(reference_groups) > 0) {
    msg <-
      paste(msg,
        "There are duplicated group numbers in the Reference group",
        sep = "\n"
      )
  }

  if (anyDuplicated(test_groups) > 0) {
    ReturnedValue$rval <- FALSE
    ReturnedValue$msg <-
      paste(ReturnedValue$msg,
        "There are duplicated group numbers in the test group",
        sep = "\n"
      )
  }
  if (is.null(model_list)) {
    msg <-
      paste(msg,
        "Model list is NULL, error in json file?",
        sep = "\n"
      )
  } else {
    for (this_model in length(model_list)) {
      if (!file.exists(model_list[this_model])) {
        msg <-
          paste(msg,
            paste0("Cannot find ", model_list[this_model]),
            sep = "\n"
          )
        next
      } else {
        control <-
          readLines(model_list[this_model], encoding = "UTF-8", warn = FALSE)

        data_line <- get_block("$DATA", control)
        if (nchar(data_line) == 0) {
          msg <- paste(msg,
            paste0("In the model file ", model_list[this_model], " DATA block not found."),
            sep = "\n"
          )
          next
        }

        data_line <-
          gsub("(\\s*\\$DATA\\s*)|(\\s*$)", "", data_line)
        any.quotes <- grep("^\"", data_line)
        if (length(any.quotes) > 0) {
          data_file <-
            sapply(
              regmatches(
                data_line,
                gregexpr('(\").*?(\")', data_line, perl = TRUE)
              ),
              function(y) {
                gsub("^\"|\"$", "", y)
              }
            )[1]
        } else {
          # find first white space
          data_file <- unlist(strsplit(data_line, " "))[1]
        }

        if (!file.exists(data_file)) {
          msg <-
            paste(msg,
              paste("Cannot find", data_file),
              sep = "\n"
            )
          next
        }

        # repeat IDs??
        Data <-
          utils::read.csv(data_file, stringsAsFactors = FALSE)
        # no ID found
        if (!"ID" %in% colnames(Data)) {
          msg <-
            append(msg, list(
              rval = FALSE,
              msg = paste("ID column not found in data set ", data_file)
            ))
        }

        IDtimesChanged <-
          which(c(FALSE, tail(Data$ID, -1) != head(Data$ID, -1)))

        if (length(IDtimesChanged) + 1 != length(unique(Data$ID))) {
          msg <-
            paste(
              msg,
              paste0(
                "There appears to be repeat IDs in ",
                data_file,
                "\n, for bootstrap sampling IDs must not repeat. "
              ),
              sep = "\n"
            )
        }

        # $EST must be on one line, need to fix this contains ;;;; Start EST?
        contains_start <-
          grepl(";;;;.*Start\\s+EST", control, ignore.case = TRUE)

        if (!any(contains_start)) {
          msg <- paste(
            msg,
            paste0(
              "The control file ",
              model_list[this_model],
              " does not contain \";;;; Start EST\", required before the $EST and the $THETA records and after $OMEGA and $SIGMA"
            ),
            sep = "\n"
          )
        }
        if (use_check_identifiable) {
          EST.line <- get_block("$EST", control)
          if (!grepl("SADDLE_RESET\\s*=\\s*1", EST.line, fixed = FALSE)) {
            msg <-
              paste(
                msg,
                paste0(
                  "Identifiability check requested, but SADDLE_RESET not set to 1 in ",
                  model_list[this_model]
                ),
                sep = "\n"
              )
          }
        }
        if (use_simulation_data) {
          if (!file.exists(simulation_data_path)) {
            msg <- paste(msg,
              paste0(
                "Cannot find simulation data set",
                simulation_data_path
              ),
              sep = "\n"
            )
          }
        }
      }
    }
    ReturnedValue <- list()
    if (nchar(msg) > 0) {
      ReturnedValue$rval <- FALSE
      ReturnedValue$msg <- paste(msg, "Exiting", "\n")
    } else {
      ReturnedValue$rval <- TRUE
      ReturnedValue$msg <- "Passes requirements check"
    }
  }
  return(ReturnedValue)
}

#' split_path
#'
#' split a path into parents
#'
#' @param  path, string, path name to be split
#' @param mustWork logical, optional, default = FALSE
#' @returns list of path parents
split_path <- function(path, mustWork = FALSE) {
  output <- c(strsplit(dirname(normalizePath(path, mustWork = FALSE)), "/|\\\\")[[1]], basename(path))
  return(output)
}
#' copy_model_files
#'
#' copy NONMEM model files from a source to a run directory
#'
#' @param  model_source folder name where subfolders (model1..modelnmodels) with model files can be found. There must be exactly one .mod file (with any stem) in the folder
#' @param run_dir Folder where models are to be run
#' @return nmodels how many models are there
#' @examples
#' \dontrun{
#' copy_model_files("c:/modelaveraging", "c:/modelaveraging/run", 4)
#' }
copy_model_files <- function(model_source, run_dir) {
  if (!file.exists(run_dir)) {
    # split run_dir, check each parent
    normpath <- normalizePath(run_dir)
    parents <- split_path(run_dir)
    nparents <- length(parents)
    cur_path <- parents[1] # should be the drive
    for (this.parent in 2:nparents) {
      cur_path <- file.path(cur_path, parents[this.parent])
      if (!file.exists(cur_path)) {
        dir.create(cur_path)
      }
    }
  }
  nmodels <- length(model_source)
  tryCatch(
    {
      for (this_model in 1:nmodels) {
        bs_dir <- file.path(run_dir, paste0("model", this_model)) # bootstrap directory
        if (dir.exists(bs_dir)) {
          unlink(bs_dir, recursive = TRUE, force = TRUE)
        }
        dir.create(bs_dir)
        modfile <- model_source[this_model]
        if (file.exists(modfile)) {
          file.copy(modfile, file.path(bs_dir, paste0("bs", this_model, ".mod")))
        } else {
          message("Cannot file model file ", modfile, ", exiting")
          return(-999)
        }
      }
      return(nmodels)
    },
    error = function(conc) {
      message("Error in copy_model_files, model # ", this_model, ", exiting")
      return(-999)
    }
  )
}


#' sample_data
#'
#' create bootstrap samples of NONMEM data set, placed in run_dir, file names = data_sampN.csv
#'
#' @param run_dir Folder where models are to be run
#' @param samp_size how many bootstrap samples
#' @param nmodels how many models are there
#' @examples
#' \dontrun{
#' sample_data("c:/modelaveraging", 100, 4)
#' }
sample_data <- function(run_dir, nmodels, samp_size) {
  con <- file(file.path(run_dir, "model1", "bs1.mod"), "r")
  suppressWarnings(control <- readLines(con, encoding = "UTF-8"))
  close(con)
  data_line <- stringr::str_trim(control[grep("\\$DATA", control)])
  data_line <- stringr::str_trim(stringr::str_replace(data_line, "\\$DATA", ""), side = "both")
  # if file name is quoted, just put out part in in quotes, otherwise get first white space
  any.quotes <- grep("^\"", data_line)
  if (length(any.quotes) > 0) {
    # find 2nd
    pos <- gregexpr(pattern = "\"", data_line)
    data_file <- substr(data_line, 1, pos[[1]][2])
    rest.of.data_line <- substr(data_line, pos[[1]][2], nchar(data_line))
  } else {
    # find first white space
    pos <- gregexpr(pattern = "\\s", data_line)
    data_file <- stringr::str_trim(substr(data_line, 1, pos[[1]][1]), side = "both")
    rest.of.data_line <- substr(data_line, pos[[1]][1], nchar(data_line))
  }
  # get datafile possibly different data files in different models? not supported at this time
  org.data <- utils::read.csv(data_file, header = TRUE, stringsAsFactors = FALSE)
  # get IDs
  cols <- colnames(org.data)
  num.data.items <- length(cols)
  # cannot have repeat IDs!!!!
  idID <- match(cols, "ID")
  IDcol <- which(idID %in% 1)
  IDs <- org.data[IDcol] %>%
    dplyr::distinct(ID)
  nsubs <- dim(IDs)[1]
  # create BS data sets
  for (this_samp in 1:samp_size) {
    who <- sample(IDs$ID, nsubs, replace = TRUE)
    this_data <- data.frame(matrix(-99, nrow = 0, ncol = num.data.items))
    colnames(this_data) <- cols
    for (this_rep in 1:nsubs) {
      # need to be sure not to have adjacent with same ID, just just number IDs sequentially
      next_sub <- org.data %>%
        dplyr::filter(ID == who[this_rep]) %>%
        dplyr::mutate(ID = this_rep)
      this_data <- rbind(this_data, next_sub)
    }
    # write data, will be in parent directory of run directory
    datafile.name <- file.path(run_dir, paste0("data_samp", this_samp, ".csv"))
    write.csv(this_data, datafile.name, quote = FALSE, row.names = FALSE)
    # replace $DATA for each control, each model

    for (this_model in 1:nmodels) {
      control_path <- file.path(run_dir, paste0("model", this_model), paste0("bs", this_model, ".mod"))
      con <- file(control_path, "r")
      suppressWarnings(control <- readLines(con, encoding = "UTF-8"))
      close(con)
      data_line <- grep("\\$DATA", control)
      # get rest of data line !!! note, will need to read entire $DATA block, maybe on more than one line!!!!
      newdata_line <- paste0("$DATA ", file.path("..", "..", paste0("data_samp", this_samp, ".csv ", rest.of.data_line)))
      control[data_line] <- newdata_line
      # remove covariance
      cov_line <- grep("\\$COV", control)
      if (length(cov_line) > 0) {
        control[cov_line] <- ";; no covariance step"
      }
      # and any tables
      control <- stringr::str_replace(control, "$TABLE", ";$TABLE")
      dir.create(file.path(run_dir, paste0("model", this_model), this_samp))

      con <- file(file.path(run_dir, paste0("model", this_model), this_samp, paste0("bsSamp", this_model, "_", this_samp, ".mod")), "w")
      writeLines(control, con)
      close(con)
    }
  }
}



#' run_one_model
#'
#' Run a single NONMEM model. This is called for either bootstrap (BS==True) or simulation,
#' The only difference in the call is the name of the resulting control, output, xml etc files
#' executable is always nonmem.exe (Windows)
#'
#' @param run_dir location of source control files, naming is bsSampModelNum_SampleNum.mod
#' @param nmfe_path path to nmfe??.bat
#' @param this_model which of the models used for model averaging is this, integer if bootstrap, NULL if monte carlo
#' @param this_samp  which sample (from bootstrap or Monte Carlo) is this, integer
#' @param BS Logical, TRUE if boostrap, FALSE if Monte Carlo
#'
#' @examples
#' \dontrun{
#'  run_one_model("c:/MBBE/rundir", "c:/nm74g64/util/nmfe74.bat", 1, 1, TRUE)
#' }
run_one_model <- function(run_dir, nmfe_path, this_model, this_samp, BS) {
  try({
    if (BS) {
      nmrundir <- file.path(run_dir, paste0("model", this_model), this_samp)
      control_file <- paste0("bsSamp", this_model, "_", this_samp, ".mod")
      output_file <- paste0("bsSamp", this_model, "_", this_samp, ".lst")
      exefile <- paste0("bsSamp", this_model, "_", this_samp, ".exe")
    } else {
      nmrundir <- file.path(run_dir, paste0("MBBEsim", this_samp))
      control_file <- paste0("MBBEsim", this_samp, ".mod")
      output_file <- paste0("MBBEsim", this_samp, ".lst")
      exefile <- paste0("MBBEsim", this_samp, ".exe")
    }

      nmoutput <- processx::run(nmfe_path, args = c(control_file, output_file), wd = nmrundir)
      delete_files(nmrundir)
  })
}
#' run_any_models
#'
#' run the bootstrap or monte carlo models/samples
#' This function can use, if available, paralllel execution with multicore or multisession
#' The plan is set in run_mbbe
#'
#' @param nmfe_path path to nmfe??.bat
#' @param run_dir Folder where models are to be run
#' @param nmodels how many models are there
#' @param samp_size how many samples are there
#' @param BS logical TRUE = is bootstrap, FALSE is Monte Carlo Simulation
#' @examples
#' \dontrun{
#' run.bootstrap("c:/nmfe744/util/nmfe74.bat", "c:/modelaveraging", 8)
#' }
run_any_models <- function(nmfe_path, run_dir, nmodels, samp_size, BS) {

  if (BS) {
    rval <- tryCatch(
      {
        total_runs <- samp_size * nmodels
        message(format(Sys.time(), digits = 0), " Starting bootstrap sample models, total number = ", total_runs)
        message("Progress bar is NONMEM bootstrap models started")
        output <- furrr::future_map(1:total_runs, function(this_num) {
          tryCatch(
            {
              this_samp <- (this_num - 1) %% samp_size + 1
              this_model <- (this_num - this_samp) / samp_size + 1
              run_one_model(run_dir, nmfe_path, this_model, this_samp, BS)

            },
            error = function(cond) {
              message("Failed to run boostrap in run_any_models, ", cond)
              return(NULL)
            }
          )
        }, .progress = TRUE, .options = furrr::furrr_options(seed = NULL))
      },
      error = function(cond) {
        message("Failed running bootstrap in run_any_models, ", cond)
        return(FALSE)
      }
    )
  } else {
    rval <- tryCatch(
      {
        total_runs <- samp_size
        message(format(Sys.time(), digits = 0), " Starting simulation models, total number = ", total_runs)
        message("Progress bar is NONMEM simulation runs started")
        output <- furrr::future_map(1:total_runs, function(this_num) {
          tryCatch(
            {
              run_one_model(run_dir, nmfe_path, NULL, this_num, FALSE)
            },
            error = function(cond) {
              message("Failed to run simulation in run_any_models, ", cond)
              return(NULL)
            }
          )
        }, .progress = TRUE, .options = furrr::furrr_options(seed = NULL))
      },
      error = function(cond) {
        message("Failed running simulation in run_any_models ", cond)
        return(NULL)
      }
    )
  }
  return(TRUE)
}

#' get_parameters
#'
#' read the xml file for each bootstrap/model, calculate BIC
#'  return list of BICs (with the BEST column) and parameters for all models/samples
#'
#' @param run_dir Folder where models are to be run
#' @param nmodels how many models are there
#' @param samp_size how many samples are there
#' @param delta_parms criteria for identifiability test
#' @param crash_value value for failed calculation
#' @param use_check_identifiable - logical, is check_identifiable to be used
#' @param passes_identifiable - array of length nmodels, count # of model that are identifable, initially can be NULL
#' @examples
#' \dontrun{
#' get_parameters("c:/modelaveraging", 4, 100, 0.1, 999999, TRUE)
#' }
get_parameters <- function(run_dir, nmodels,
                           samp_size, delta_parms,
                           crash_value, use_check_identifiable,
                           passes_identifiable) {
  BICS <- data.frame(matrix(crash_value, nrow = samp_size, ncol = nmodels + 3))
  colnames(BICS) <- c(paste0("Model", seq(1:nmodels)), "Best", "Max_Delta_parm", "Max_Delta")
  BICS$Best <- NA
  # BIC=k*ln(n) -2LL where k is is the number of estimated parameters and n is the number of data
  nparms <- data.frame(ntheta = rep(NA, nmodels))

  # only need parameters for best model, but need to read xml file to get the number of observations and parameters, so may as well get
  # THETA etc now for all
  nfailed_ident <- 0 # number of samples in this model that fail the identifiability test
  # do by sample first only save parameters for the best model
  selected_parameters <- vector("list", samp_size)
  tryCatch(
    {
      for (this_samp in 1:samp_size) {
        num_successful <- 1 # only save parameters if model finished, doesn't need to converge, just finish, this is the number of successful samples for this model
        parameters_this_sample <- vector("list", nmodels)
        this_model <- 1
        all_identifiables <- rep(TRUE, nmodels)
        for (this_model in 1:nmodels) {
          theta <- NA
          xml_file <- file.path(run_dir, paste0("model", this_model), this_samp, paste0("bsSamp", this_model, "_", this_samp, ".xml"))
          identifiable_ok <- c(passes = FALSE, max_delta = -999, max_delta_parm = -999)
          rval <- tryCatch(
            {
              count <- 0
              boolFalse <- FALSE
              # wait for files to close?? otherwise get Error in file(file, "rt") : cannot open the connection
              while (boolFalse == FALSE & count < 10) {
                data <- NULL
                count <- count + 1
                tryCatch({
                  data <- xml2::read_xml(xml_file, encoding = "ASCII")
                  boolFalse <- TRUE
                }, error = function(e) {
                  Sys.sleep(0.1)
                }, finally = {}) # fails will throw error below in outer tryCatch
              }
              if (count >= 10) {
                message("Failed to read xml file for sample ", this_samp, ", model = ", this_model, ", xml file = ", xml_file)
                BICS[this_samp, this_model] <- crash_value
                identifiable_ok <- c(passes = FALSE, max_delta = -999, max_delta_parm = -999)
                parameters_this_sample[[this_model]] <- NA
                tryCatch(
                  {
                    lstFile <- file.path(run_dir, paste0("model", this_model), this_samp, paste0("bsSamp", this_model, "_", this_samp, ".lst"))
                    cat(BICS[this_samp, this_model], file = lstFile, append = TRUE)
                  },
                  error = function(err) {
                    con <- file(lstFile, "a")
                    message("Cannot write to output file for sample, ", this_samp, ", model, ", this_model)
                  }
                )
              } else {
                problem_node <- xml2::xml_find_all(data, "//nm:problem_information")
                contents <- xml2::xml_contents(problem_node)
                text <- xml2::xml_text(contents)
                text <- as.list(unlist(strsplit(text, "\n")))
                nobsline <- grep("TOT. NO. OF OBS RECS:", text)
                nobs <- as.integer(stringr::str_replace(text[nobsline], "TOT. NO. OF OBS RECS:", ""))
                info_node <- xml2::xml_find_all(data, "//nm:problem_options")
                info_contents <- xml2::xml_attrs(info_node)
                ntheta <- as.numeric(info_contents[[1]]["nthetat"])
                if (this_samp == 1) {
                  # have to do this here, don't know how many thetas until here
                  nparms$ntheta[this_model] <- ntheta
                }
                parameters_this_model <- data.frame(matrix(NA, ncol = ntheta))
                colnames(parameters_this_model) <- paste0("THETA", seq(1:ntheta))
                parameters_this_sample[[this_model]] <- parameters_this_model
                # and OFV
                estim.node <- xml2::xml_find_all(data, "//nm:estimation")
                estim_contents <- xml2::xml_attrs(estim.node)
                # nm:final_objective_function
                OFV.node <- xml2::xml_find_all(data, "//nm:final_objective_function")
                OFV_contents <- xml2::xml_contents(OFV.node)
                OFV <- as.numeric(xml2::xml_text(OFV_contents))
                if (length(OFV) > 0) {
                  BICS[this_samp, this_model] <- ntheta * log(nobs) + OFV
                } else {
                  BICS[this_samp, this_model] <- crash_value
                }
                theta_node <- xml2::xml_find_all(data, "//nm:theta")
                theta_children <- xml2::xml_children(theta_node)
                theta_contents <- xml2::xml_contents(theta_node)
                theta <- as.numeric(xml2::xml_text(theta_contents))
                # length of theta will be 1 if a crash
                if (length(theta) > 1) {
                  parameters_this_model <- theta
                  num_successful <- num_successful + 1
                  parameters_this_sample[[this_model]] <- parameters_this_model
                } else {
                  parameters_this_sample[[this_model]] <- NA
                }

                # if using saddle_reset, test for identifiability
                if (!use_check_identifiable) {
                  identifiable_ok["passes"] <- TRUE
                  passes_identifiable[this_model] <- passes_identifiable[this_model] + 1
                } else {
                  # run_dir,this_model,this_sample,delta_parms
                  identifiable_ok <- check_identifiable(run_dir, this_model, this_samp, delta_parms, ntheta)
                  all_identifiables[this_model] <- identifiable_ok["passes"]
                  if (!identifiable_ok$passes) {
                    # first just get the number of parameters and observations for this model
                    BICS[this_samp, this_model] <- crash_value
                  } else {
                    passes_identifiable[this_model] <- passes_identifiable[this_model] + 1
                  }
                }
              }
            },
            error = function(e) {
              message("error in get_parameters, read data, inner loop, error  = ", cond, ", sample ", this_samp, ", model = ", this_model)
              identifiable_ok <- c(passes = FALSE, max_delta = -999, max_delta_parm = -999)
              all_identifiables[this_model] <- identifiable_ok["passes"]
              BICS[this_samp, this_model] <- crash_value
              parameters_this_sample[[this_model]] <- NA
            }
          )
          tryCatch(
            {
              lstFile <- file.path(
                run_dir, paste0("model", this_model),
                this_samp, paste0("bsSamp", this_model, "_", this_samp, ".lst")
              )
              cat(
                paste0(
                  "BIC = ", round(BICS[this_samp, this_model], 3),
                  ", identifiable = ", identifiable_ok["passes"],
                  ", max_delta = ", round(identifiable_ok$max_delta, 5)
                ),
                file = lstFile, append = TRUE
              )
              cat(paste0(", max_delta_parm = ", identifiable_ok$max_delta_parm),
                file = lstFile, append = TRUE
              )
            },
            error = function(err) {
              message("Failed to append parameters to ", lstFile, " after failing in get_parameters")
            }
          )
        } # end of model loop
        #  }  # end of sample loop

        # and select best model, NA if all crashed
        if (all(BICS[this_samp, 1:nmodels] == crash_value)) {
          BICS$Best[this_samp] <- NA
          warning("No models successful for model ", this_samp, " identifiability = ", toString(all_identifiables))
          selected_parameters[[this_samp]] <- NA
          BICS$Max_Delta[this_samp] <- 999999
          BICS$Max_Delta_parm[this_samp] <- NA
        } else {
          best <- which.min(BICS[this_samp, 1:nmodels])
          BICS$Best[this_samp] <- best
          selected_parameters[[this_samp]] <- unlist(parameters_this_sample[[best]])
          # check in nparms = length of selected_parameters
          message("\nSample ", this_samp, " best model is ", best, " with parameters:\n")

          for (i in 1:length(selected_parameters[[this_samp]])){
            message("parameter[", i, "] = ", round(selected_parameters[[this_samp]][i],6))
          }
          if(nparms$ntheta[best]!= length(selected_parameters[[this_samp]])){
            stop("Length of parameters for sample ", this_samp, " doesn't equal the number in the xml file for model ", best, ", exiting")
          }
          BICS$Max_Delta[this_samp] <- identifiable_ok["max_delta"]
          BICS$Max_Delta_parm[this_samp] <- identifiable_ok["max_delta_parm"]
        }
      }
    },
    error = function(cond) {
      # everything crashes??
      message("Unknown error in get_parameters, read data, sample ", this_samp, " error  = ", cond)
      all_identifiables[this_model] <- identifiable_ok["passes"]
      BICS[this_samp, this_model] <- crash_value
    }
  )
  # write out results,
  tryCatch(
    {
      flatBICs <- apply(BICS, 2, as.character)
      write.csv(flatBICs, file.path(run_dir, "BICS.csv"), quote = FALSE)
      message("BICS = ")
      message(paste0(capture.output(BICS), collapse = "\n"))
    },
    error = function(err) {
      message("Failed to write out BICS to ", file.path(run_dir, "BICS.csv"), " in get_parameters")
    }
  )
  # write out parameters, can't do write.csv, as rows have different number of parameters

  tryCatch(
    {
      conn <- file(file.path(run_dir, "Parameters.csv"), open = "w")

      writeLines(text = c("Sample","Best Model for Sample","Number of Parmeters","Selected Parameters"), con = conn, sep = "\n")
      newlist <- lapply(seq_len(length(selected_parameters)), function(i) {
        temp <- lapply(seq_len(length(selected_parameters[[i]])), function(j) {
          temp <- c(selected_parameters[[i]][[j]])
        })

        writeLines(text = paste(i,BICS$Best[i], length(selected_parameters[[i]]),
                                paste(temp, collapse = ","), sep = ","), con = conn, sep = "\n")
      })
      close(conn)
    },
    error = function(err) {
      message("Failed to write to ", file.path(run_dir, "Parameters.csv"), " in get_parameters")
    }
  )
  rval <- list(BICS = BICS, parameters = selected_parameters, passes_identifiable = passes_identifiable)
  return(rval)
}

#' get_base_model
#'
#' for each model get the control file used for the bootstrap
#' return a list of the models
#' @param run_dir directory where models are run
#' @param nmodels how many models are there
#'
#' @examples
#' \dontrun{
#' get_base_model("c:/mbbe/rundir", 4)
#' }
get_base_model <- function(run_dir, nmodels) {
  # need error trapping for no ;;;;;.*Start EST
  base_models <- vector(mode = "list", length = 0) # all but $THETA, $SIM, $TABLE

  for (this_model in 1:nmodels) {
    tryCatch(
      {
        con <- file(file.path(run_dir, paste0("model", this_model), paste0("bs", this_model, ".mod")), "r")
        suppressWarnings(control <- readLines(con, encoding = "UTF-8"))
        control <- control[1:grep(";;;;.*Start EST", control)]
        close(con)
        base_models <- append(base_models, list(control))
      },
      error = function(e) {
        stop("Cannot find ", file.path(paste0(run_dir, "model", this_model), paste0("bs", this_model, ".mod")), "exiting")
      }
    )
  }
  return(base_models)
}

#' Edits the best based model for each sample, replaces the original parameters with the bootstrap parameters
#' and adds the $SIM and $TABLE
#'
#' @param run_dir Folder where models are to be run
#' @param parms list that include BICs and parameters for each sample/model
#' @param parms list that include BICs and parameters for each sample/model
#' @param base_models list of the text in each model used for the bootstrap
#' @param samp_size how many samples are there
#' @param use_simulation_data logical - TRUE is a file other than the bootstrap data file is to be used for simulation
#' @param simulation_data_path if use_simulation_data==TRUE, the path to the file tto be used for simulation
#' @examples
#' \dontrun{
#' write_sim_controls("c:/modelaveraging", parms, base_models, 100)
#' }
write_sim_controls <- function(run_dir, parms, base_models, samp_size, use_simulation_data, simulation_data_path = NULL) {
  final_models <- list()
  if (!file.exists(simulation_data_path)) {
    stop(paste("Cannot find ", simulation_data_path, " exiting"))
  } else {
    nmodels <- length(base_models)
   # mbfinal_models <- vector(mode = "list", length = samp_size)
    #model.indices <- rep(0, nmodels) # which model and parameter set to use when this model is selected, roll over if not enough samples,
    current_runable_samp <- 0 # some models may not be runable (all crashed), if so, start over with model, and different seed in NONMEM
    for (this_samp in 1:samp_size) {
      current_runable_samp <- current_runable_samp + 1
      count <- 0
      while (is.na(parms$BICS$Best[current_runable_samp]) & count < samp_size * 4) { # only go through samples 4 times? don't keep running same samples over?
        current_runable_samp <- current_runable_samp + 1
        count <- count + 1
        if (current_runable_samp > samp_size) {current_runable_samp <- 1}
        # use BS parameters, when you run out (as some BS samples fail), just start over, so need different random seed in $SIM model
      }
      if (count >= (samp_size * 4)) {
        stop("Not enough samples with successful outcome for simulation")
      }
      # current_runable_samp is sample #, from 1 to samp_size
      # which.model is the best model in that sample, from 1 to nmodels
      which.model <- parms$BICS$Best[current_runable_samp]
      # is first index, e.g., parameters[[2]]$THETA1[1] is THETA[1] for model 2, first sample
      # model.indices[which.model] <- model.indices[which.model] + 1
      full_control <- base_models[which.model][[1]]
      # need to get sequence in here, calculated in $PK
      seed <- round(runif(1, 0, 10000), 0)
      full_control <- c(full_control, paste("$SIM ONLYSIM (", seed, ")"), "$THETA")
      ntheta <- length(parms$parameters[[which.model]])

      if (use_simulation_data) {
        # get $DATA
        start_line <- grep("^\\$DATA", full_control)
        next_line <- grep("^\\$", full_control[start_line + 1:length(full_control)])[1]
        if (is.empty(next_line)) {
          # $DATA is last
          last_line <- length(full_control[[1]])
        } else {
          last_line <- start_line + next_line
        }
        # replace with simulation data set
        full_control <- c(full_control[1:(start_line - 1)], paste("$DATA", simulation_data_path, "IGNORE=@"), full_control[(last_line):length(full_control)])
      }
      for (this_parm in parms$parameters[[which.model]]) {
        full_control <- c(full_control, paste0(this_parm, "  ;; THETA(", this_parm, ")"))
      }
      full_control <- c(full_control, "$TABLE ID TIME GROUP OCC SEQ DV EVID NOPRINT NOAPPEND FILE=OUT.DAT ONEHEADER")
      full_control <- c(full_control, paste0(";;Sample #", this_samp, "; Selected sample = ",current_runable_samp, "\n"))
      full_control <- c(full_control, paste0(";;Source = model #", which.model, ", numparms = ", ntheta, "\n"))

      sim_dir <- file.path(run_dir, paste0("MBBEsim", this_samp))
      if (dir.exists(sim_dir)) {
        unlink(sim_dir, recursive = TRUE, force = TRUE)
      }
      if (dir.exists(sim_dir)) {
        unlink(sim_dir, recursive = TRUE, force = TRUE)
      }
      dir.create(sim_dir)
      control_file <- file.path(sim_dir, paste0("MBBEsim", this_samp, ".mod"))
      count <- 0
      while (!file.exists(control_file) & count < 10) {
        con <- file(control_file, "w")
        Sys.sleep(0.1)
        writeLines(unlist(full_control), con)
        count <- count + 1
        close(con)
      }
      if (!file.exists((control_file))) {
        message("Cannot write to ", control)
      }
      final_models <- append(final_models, list(full_control))
    }
    return(final_models)
  }
}
#'

#' calc_NCA
#'
#' perform NCA (Cmax, AUCin,AUClst, from 0 to end.time)
#' @param run_dir character;  folder path of run directory
#' @param ngroups numeric; integer value specifying how many groups, e.g., 4 for ABBA
#' @param reference_groups numeric vector; specifying which groups are reference
#' @param test_groups numeric vector; specifying which groups are test
#' @param NCA_end_time numeric; end time for AUClast and AUCinf
#' @param samp_size numeric; integer value specifying sample size
#' @examples
#' \dontrun{
#' run_dir <- "c:/Workspace/mbbe"
#' ngroups <- 4
#' reference_groups <- c(1,2)
#' samp_size <- 6
#' test_groups <- c(3,4)
#' NCA_end_time <- 7
#' calc_NCA(run_dir, ngroups, reference_groups, test_groups, NCA_end_time, samp_size)
#' }
#' @details
#' calc_NCA calls getNCA for each sample in 1:samp_size
#' Currently not parallelized
#'
#'
#'@export
calc_NCA <-
  function(run_dir,
           ngroups,
           reference_groups,
           test_groups,
           NCA_end_time,
           samp_size) {
    # still need to parallelize
    tryCatch(
      {
        for (this_samp in 1:samp_size) {
          getNCA(
            run_dir,
            this_samp,
            ngroups,
            reference_groups,
            test_groups,
            NCA_end_time
          )
        }
      },
      error = function(cond) {
        message("Failed in calc_NCA, sample # ", this_samp, ", error ", cond)
        return(NULL)
      }
    )
    return()
  }

#' check_identifability
#'
#' Identifability is defined by Aoki (https://www.page-meeting.org/default.asp?abstract=5951)
#' Briefly, a saddle_reset is required for the NONMEM minimzation. The pre-reset parameters are compared to the
#' post-reset parameters. If the fractional difference between any pre and post reset parameter exceeds delta_parms
#' the model fails identifiability
#' Algorithm: reads .lst file, find out where saddle reset occurs (doesn't seem to be recorded anywhere else)
#' then get parameter before reset from .ext file (Pre.parms)
#' compare with final parameters (Post.parms), see if any differ by delta_parms
#'
#' @param run_dir home directory
#' @param this_model integer
#' @param this_sample integer
#' @param delta_parms absolute difference in parameters criteria for failing identifability
#' @param nparms number of parameters in the .ext file (can be less than the total)
#' @return
#' returns a list consisting of:
#' passes - logical
#' max_delta - maximum difference seen in pre and post reset parameters
#' max_delta_parm which parameter in the .ext file have the largest difference
#' @examples
#' \dontrun{
#' check_identifiable("c:/runmodels", 1, 1, 0.1)
#' }
#'
check_identifiable <- function(run_dir, this_model, this_sample, delta_parms, nparms) {
  lstfile <- file.path(run_dir, paste0("model", this_model), this_sample, paste0("bsSamp", this_model, "_", this_sample, ".lst"))
  extfile <- file.path(run_dir, paste0("model", this_model), this_sample, paste0("bsSamp", this_model, "_", this_sample, ".ext"))
  max_deltap <- -999
  max_delta_parmp <- -999
  if (!file.exists(lstfile) | !file.exists(extfile)) {
    message("Cannot find ", listfile, " or ", extfile, " for determining identifiability")
    return(c(passes = FALSE, max_delta = -999, max_delta_parm = -999))
  } else {
    tryCatch(
      {
        con <- file(lstfile, "r")
        suppressWarnings(output <- readLines(con, encoding = "UTF-8"))
        close(con)
        reset.line <- grep("^0SADDLE POINT RESET", output)

        # get previous '0ITERATION NO.: '
        first_output <- output[1:reset.line]
        first_output <- first_output[grep("^0ITERATION NO.:", first_output)]
        first_output <- first_output[length(first_output)]
        reset_iteration <- as.integer(substr(first_output, 16, 24))
        last_output <- output[reset.line:length(output)]
        last_output <- last_output[grep("^0ITERATION NO.:", last_output)]
        last_output <- last_output[1]
        last.iteration <- as.integer(substr(last_output, 16, 24))
        # read parameters from .ext
        ext <- read.table(extfile, header = TRUE, skip = 1)
        Pre.parms <- ext %>%
          dplyr::filter(ITERATION == reset_iteration)
        Post.parms <- ext %>%
          dplyr::filter(ITERATION == last.iteration)
        # nparms <- length(Pre.parms)  # includes the first column - 'ITERATION'
        Passes_identifiability <- TRUE

        for (this_parm in 2:nparms) {
          if (Post.parms[this_parm] != 0) {
            difference <- abs((Pre.parms[this_parm] - Post.parms[this_parm]) / Post.parms[this_parm])
            if (difference > delta_parms) {
              Passes_identifiability <- FALSE
            }
            if (difference > max_deltap) {
              max_deltap <- difference
              max_delta_parmp <- this_parm - 1
            }
          }
        }
        rval <- c(passes = Passes_identifiability, max_delta = max_deltap, max_delta_parm = max_delta_parmp)
        names(rval) <- c("passes", "max_delta", "max_delta_parm")
        return(rval)
      },
      error = function(cond) {
        rval <- c(passes = Passes_identifiability, max_delta = max_deltap, max_delta_parm = max_delta_parmp)
        names(rval) <- c("passes", "max_delta", "max_delta_parm")
        return(rval)
      }
    )
  }
}

#' getNCA
#'
#' Calculated NCA parameters (Cmax, AUCinf, AUClast) by group from NONMEM output file
#' NONMEM output file MUST be in the folder MBBEsimN where N is the simulation number
#' $TABLE must include the ONEHEADER option
#' The simulation must include an observation at T=
#' The $TABLE file output by the simulation must include:
#' ID
#' GROUP
#' SEQ
#' OCC
#' Note that for MBBE, the simulation control file is written by the write_sim function and will include:
#'  "$TABLE ID TIME GROUP OCC SEQ DV EVID NOPRINT NOAPPEND FILE=OUT.DAT ONEHEADER")
#' reads $TABLE output from simulation (file name out.dat)
#' and do NCA (Cmax, AUCin,AUClst, from 0 to end.time)
#' if check.identifiability, do that with delta_parms as criteria
#' @param run_dir home directory
#' @param this_sample integer
#' @param NumGroups integer how many groups, e.g., 4 for ABBA
#' @param reference_groups list of arrays, which groups are reference
#' @param test_groups list of arrays, which groups are  test
#' @param NCA_end_time end time for AUClast and AUCinf
#' @examples
#' \dontrun{
#' getNCA("c:/runmodels", 1, 4, c(1, 2), c(3, 4), 72, 0.1)
#' }
getNCA <- function(run_dir, this_sample, NumGroups, reference_groups, test_groups, NCA_end_time) {
  output_file <- file.path(run_dir, paste0("MBBEsim", this_sample), paste0("NCAresults", this_sample, ".csv"))
  tryCatch(
    {
      NMoutFile <- file.path(run_dir, paste0("MBBEsim", this_sample), "out.dat")
      if (file.exists(NMoutFile)) {
        group_NCA_results <- data.frame(
          ID = as.integer(), treatment = as.integer(), period = as.integer(), sequence = as.integer(),
          Cmax = as.numeric(), AUCinf = as.numeric(), AUClast = as.numeric()
        )
        All_NCA_results <- data.frame(
          ID = as.integer(), treatment = as.integer(), period = as.integer(), sequence = as.integer(),
          Cmax = as.numeric(), AUCinf = as.numeric(), AUClast = as.numeric()
        )
        data <- read.table(NMoutFile, skip = 1, header = TRUE)
        data <- data %>%
          dplyr::filter(EVID == 0) %>%
          dplyr::filter(DV > 0)

        if (file.exists(output_file)) file.remove(output_file)
        for (this_group in 1:NumGroups) {
          group_data <- data %>%
            dplyr::filter(GROUP == this_group)
          # keep period for each subject, for this group
          period_seq <- group_data %>%
            dplyr::group_by(ID) %>%
            dplyr::distinct(ID, .keep_all = TRUE) %>%
            dplyr::select(ID, OCC, SEQ) %>%
            dplyr::arrange(ID)

          # insert conc=0 at time = 0, but remove if duplicated??
          zero_time <- group_data %>%
            dplyr::distinct(ID, .keep_all = TRUE)
          zero_time$TIME <- 0
          zero_time$DV <- 0
          group_data <- rbind(group_data, zero_time) %>%
            dplyr::arrange(ID, TIME)
          conc_obj <- PKNCA::PKNCAconc(
            group_data,
            DV ~ TIME | ID
          )
          data_obj <- PKNCA::PKNCAdata(
            data.conc = conc_obj,
            intervals = data.frame(
              start = 0,
              end = NCA_end_time,
              aucinf.obs = TRUE,
              auclast = TRUE,
              cmax = TRUE
            )
          )
          suppressMessages(
            results_obj <- PKNCA::pk.nca(data_obj, verbose = FALSE)$result
          )
          AUCinf <- results_obj %>%
            dplyr::filter(PPTESTCD == "aucinf.obs") %>%
            dplyr::select(ID, PPORRES)
          AUClast <- results_obj %>%
            dplyr::filter(PPTESTCD == "auclast") %>%
            dplyr::select(ID, PPORRES)
          CMAX <- results_obj %>%
            dplyr::filter(PPTESTCD == "cmax") %>%
            dplyr::select(ID, PPORRES)
          if (this_group %in% reference_groups) {
            treatment <- "Reference"
          } else {
            treatment <- "Test"
          }
          group_NCA_results <- data.frame(
            ID = AUCinf$ID,
            treatment = treatment,
            period = period_seq$OCC,
            sequence = period_seq$SEQ,
            Cmax = CMAX$PPORRES,
            AUCinf = AUCinf$PPORRES,
            AUClast = AUClast$PPORRES
          )
          All_NCA_results <- rbind(
            All_NCA_results,
            group_NCA_results
          )
        }
        # wait for file to be written??, this doesn't seem to help


        count <- 0
        while (file.exists(output_file) & count < 50) {
          file.remove(output_file)
          count <- count + 1
          Sys.sleep(0.1)
        }
        count <- 0
        while (!file.exists(output_file) & count < 50) {
          write.csv(All_NCA_results,
            file = output_file,
            quote = FALSE,
            row.names = FALSE
          )
          count <- count + 1
          Sys.sleep(0.1)
        }
        if (!file.exists(output_file)) {
          message("Cannot write to ", output_file)
          return(FALSE)
        }
        return(TRUE)
      } else {
        return(FALSE)
      }
    },
    error = function(cond) {
      message("Error in NCA calculation, ", this_sample, "\n")
      message(cond)
      return(FALSE)
    }
  )
}

#' Calculate power
#'
#' Calculated power by doing EMA standards statistics on each Monte Carlo simulation, then counting the number that pass BE
#' @param run_dir character, run directory
#' @param samp_size integer, how many samples
#' @param alpha numeric, range, >0 and < 1 alpha error rate
#' @param model_averaging_by character "subject" or "study"
#' @param NTID logical, is this narrow therapeutic index drug
#' @export
#' @return a list of Cmax_result, AUCinf_result, AUClast_result power (0-1)
#' @details
#' The simulation is done by study (e.g., one model per study) this will result in model averaging
#' by study. If model_averaging_by == "subject", the by-study data are recombined, randomly chosing one subject
#' (without replacement) from all studies to reassemble each study data set
#' The function loops over each sample, read the NCAresultsN where N is the sample number
#' Then calculates whether that sample passes or fails BE testing
#' @examples
#' \dontrun{
#' nca_results <- system.file(package = "mbbe", "nca_results")
#' calc_power("c:/MBBErundir", 100, 0.05, "study", FALSE)
#' }

calc_power <- function(run_dir, samp_size, alpha, model_averaging_by, NTID) {

  BEsuccess <- data.frame(
    Cmax.success = as.logical(),
    AUCinf.success = as.logical(),
    AUClast.success = as.logical()
  )
  message(format(Sys.time(), digits = 0), " Starting statistics for NCA parameters, simulations 1-", samp_size)
  all_results <- data.frame(
    SampleNum = as.integer(),
    Cmax_Ratio = as.numeric(),
    Cmax_lower_CL = as.numeric(),
    Cmax_upper_CL = as.numeric(),
    Cmax_swR = as.numeric(),
    Cmax_pe = as.numeric(),
    Cmax_critbound = as.numeric(),
    Cmax_Assessment = as.numeric(),
    Cmax_BE = as.logical(),
    AUCinf_Ratio = as.numeric(),
    AUCinf_lower_CL = as.numeric(),
    AUCinf_upper_CL = as.numeric(),
    AUCinf_swR = as.numeric(),
    AUCinf_pe = as.numeric(),
    AUCinf_critbound = as.numeric(),
    AUCinf_Assessment = as.numeric(),
    AUCinf_BE = as.logical(),
    AUClast_Ratio = as.numeric(),
    AUClast_lower_CL = as.numeric(),
    AUClast_upper_CL = as.numeric(),
    AUClast_swR = as.numeric(),
    AUClast_pe = as.numeric(),
    AUClast_critbound = as.numeric(),
    AUClast_Assessment = as.numeric(),
    AUClast_BE = as.logical()
  )

  # read in all NCAresults file
  # theses will one one model per study
  # once all are collected reassmble into on model per subject in each study
  all_nca_data <- data.frame(
    sample = as.integer(), ID = as.integer(), treatment = as.integer(), period = as.integer(), sequence = as.integer(), Cmax = as.numeric(),
    AUCinf = as.numeric(), AUClast = as.numeric()
  )
  nsubs <- NA
  for (this_samp in 1:samp_size) {
    this_ncafile <- file.path(run_dir, paste0("MBBEsim", this_samp), paste0("NCAresults", this_samp, ".csv"))
    # wait for files to close?? otherwise get Error in file(file, "rt") : cannot open the connection
    data <- NULL
    boolFalse <- FALSE
    count <- 0
    while (boolFalse == FALSE & count < 20) {
      data <- NULL
      count <- count + 1
      tryCatch({
        if (file.exists(this_ncafile)) {}
        this_data <- utils::read.csv(this_ncafile, header = TRUE)
        boolFalse <- TRUE
        if (is.na(nsubs)) {
          IDs <- this_data %>% dplyr::distinct(ID)
          nsubs <- dim(IDs)[1]
        }
      }, error = function(e) {
        Sys.sleep(0.1)
      }, finally = {}) # fails will throw error below
    }
    if(exists("this_data")){
      this_data$sample <- this_samp
    }else{
      stop("Cannot fine ",this_ncafile, " NCA output file")
    }

    this_data <- subset(this_data, select = c(sample, ID:AUClast))
    all_nca_data <- rbind(all_nca_data, this_data)
  }

  if (model_averaging_by == "subject") {
    # original, doesn't really need to be study numbers
    sequences <- matrix(NA, nrow = nsubs, ncol = samp_size)
    colnames(sequences) <- paste0("Study", 1:samp_size)
    for (this_sub in 1:nsubs) {
      sequences[this_sub, ] <- sample(1:samp_size, samp_size, replace = FALSE)
    }
    # reassemble all_nca_data
    #
    new_all_nca_data <- all_nca_data[0, ]
    # and fill in with new sequence
    for (this_samp in 1:samp_size) {
      for (this_sub in 1:nsubs) {
        new_sample <- all_nca_data %>% dplyr::filter(ID == this_sub, sample == sequences[this_sub, this_samp])
        new_sample$sample <- this_samp
        new_all_nca_data <- rbind(new_all_nca_data, new_sample)
      }
    }
    all_nca_data <- new_all_nca_data
    outfile <- file.path(run_dir, "resampledNCA.csv")
    message("Resampled NCA data for subject level model averaging written to ", outfile)
    write.csv(all_nca_data, outfile)
  }

  ## and do stats

  message(format(Sys.time(), digits = 0), " Done reading NCA results, doing stats")
  if (!is.null(all_nca_data)) {
    Cmax_result <- data.frame(
      MetricColumn = "Cmax", Ratio = -999, lower.CL = -999, upper.CL = -999,
      swR = -999, pe = -999, critbound = -999, VarianceCriterion = -999,
      Assessment = -999, BE = -999
    )
    tryCatch(
      {
        Cmax_result <- get_BEQDF(all_nca_data %>% dplyr::filter(!is.na(Cmax)),
          MetricColumn = "Cmax",
          SubjectColumn = "ID", TreatmentColumn = "treatment", SequenceColumn = "sequence",
          PeriodColumn = "period", RefValue = "Reference", alpha = alpha, PartialReplicate = FALSE, NTID
        )
        Cmax_result$SampleNum <- this_samp
        # reorder, but sampleNum first
        Cmax_result <- subset(Cmax_result, select = c(SampleNum, MetricColumn:BE))
      },
      error = function(e) {
        Cmax_result <- data.frame(
          SampleNum = this_samp, MetricColumn = "Cmax", Ratio = -999, lower.CL = -999, upper.CL = -999,
          swR = -999, pe = -999, critbound = -999, Assessment = -999, BE = -999
        )
      }
    )
    colnames(Cmax_result) <- c(
      "SampleNum", "Cmax_MetricColumn", "Cmax_Ratio", "Cmax_lower_CL", "Cmax_upper_CL", "Cmax_swR", "Cmax_pe",
      "Cmax_critbound", "Cmax_VarianceCriterion", "Cmax_Assessment", "Cmax_BE"
    )


    AUClast_result <- data.frame(
      MetricColumn = "AUClast", Ratio = -999, lower.CL = -999, upper.CL = -999,
      swR = -999, pe = -999, critbound = -999, VarianceCriterion = -999,
      Assessment = -999, BE = -999
    )
    tryCatch(
      {
        AUClast_result <- get_BEQDF(all_nca_data %>% dplyr::filter(!is.na(AUClast)),
          MetricColumn = "AUClast",
          SubjectColumn = "ID", TreatmentColumn = "treatment", SequenceColumn = "sequence",
          PeriodColumn = "period", RefValue = "Reference", alpha = alpha, PartialReplicate = FALSE, NTID
        )
      },
      error = function(e) {
        AUClast_result <- data.frame(
          MetricColumn = "AUClast", Ratio = -999, lower.CL = -999, upper.CL = -999,
          swR = -999, pe = -999, critbound = -999, VarianceCriterion = -999, Assessment = -999, BE = -999
        )
      }
    )
    colnames(AUClast_result) <- c(
      "AUClast_MetricColumn", "AUClast_Ratio", "AUClast_lower_CL", "AUClast_upper_CL", "AUClast_swR", "AUClast_pe",
      "AUClast_critbound", "AUClast_VarianceCriterion", "AUClast_Assessment", "AUClast_BE"
    )


    AUCinf_result <- data.frame(
      MetricColumn = "AUCinf", Ratio = -999, lower.CL = -999, upper.CL = -999,
      swR = -999, pe = -999, critbound = -999, VarianceCriterion = -999,
      Assessment = -999, BE = -999
    )
    tryCatch(
      {
        AUCinf_result <- get_BEQDF(all_nca_data %>% dplyr::filter(!is.na(AUCinf)),
          MetricColumn = "AUCinf",
          SubjectColumn = "ID", TreatmentColumn = "treatment", SequenceColumn = "sequence",
          PeriodColumn = "period", RefValue = "Reference", alpha = alpha, PartialReplicate = FALSE, NTID
        )
      },
      error = function(e) {
        AUCinf_result <- data.frame(
          MetricColumn = "AUCinf", Ratio = -999, lower.CL = -999, upper.CL = -999,
          swR = -999, pe = -999, critbound = -999, VarianceCriterion = -999,
          Assessment = -999, BE = -999
        )
      }
    )

    colnames(AUCinf_result) <- c(
      "AUCinf_MetricColumn", "AUCinf_Ratio", "AUCinf_lower_CL", "AUCinf_upper_CL", "AUCinf_swR", "AUCinf_pe",
      "AUCinf_critbound", "AUCinf_VarianceCriterion", "AUCinf_Assessment", "AUCinf_BE"
    )
    all_results <- rbind(all_results, c(Cmax_result, AUCinf_result, AUClast_result))
  } else {
    # no out.dat
  }
  return(all_results)
}
#' make_NCA_plots
#'
#'generated histograms of AUClast, AUCinf and Cmax for the BE samples
#' @param BICS - data objects for Bayesian information Criteria
#' @param run_dir
#' @param samp_size
#' @param nmodels
#' @param reference_groups
#' @param test_groups
#' @param saveplots
#'
make_NCA_plots <- function(BICS, run_dir, samp_size, nmodels, reference_groups, test_groups, saveplots = FALSE) {
  rval <- tryCatch(
    {
      BICS <- BICS %>%
        dplyr::mutate(Samp_num = dplyr::row_number())
      this_NCAs <- NULL
      all_NCAs <- data.frame(
        ID = as.integer(), treatment = as.character(), period = as.integer(), sequence = as.integer(),
        Cmax = as.numeric(), AUCinf = as.numeric(), AUClast = as.numeric(), model = as.integer()
      )
      this_model <- 1
      for (this_model in 1:nmodels) {
        all_NCAs_this_model <- data.frame(
          ID = as.integer(), treatment = as.character(), period = as.integer(), sequence = as.integer(),
          Cmax = as.numeric(), AUCinf = as.numeric(), AUClast = as.numeric(), model = as.integer()
        )
        # gather all Cmax etc from models labeled by model superimpose distributions
        which_models <- BICS$Best
        # compile list of NCA with best
        this_samp <- which_models[1]
        for (this_samp in which_models) {
          filename <- file.path(run_dir, paste0("MBBEsim", this_samp), paste0("NCAresults", this_samp, ".csv"))
          if (file.exists(filename)) {
            this_NCAs <- utils::read.csv(filename)
            this_NCAs$model <- this_model
            all_NCAs_this_model <- rbind(all_NCAs_this_model, this_NCAs)
          }
        }
        if (dim(all_NCAs_this_model)[1] > 0) {
          all_NCAs <- rbind(all_NCAs, all_NCAs_this_model)
          Cmax_plot <- ggplot2::ggplot(all_NCAs_this_model, ggplot2::aes(x = Cmax), ) +
            ggplot2::geom_histogram(ggplot2::aes(color = treatment, fill = treatment), position = "identity", alpha = 0.4, bins = 30) +
            ggplot2::ggtitle(paste("Cmax, model =", this_model))

          ggplot2::ggsave(file.path(run_dir, paste0("Model_", this_model, "_Cmax_histogram_by_treatment.jpeg")), Cmax_plot, device = "jpeg", width = 9, height = 6)
          AUCinf_plot <- ggplot2::ggplot(all_NCAs_this_model, ggplot2::aes(x = AUCinf), ) +
            ggplot2::geom_histogram(ggplot2::aes(color = treatment, fill = treatment), position = "identity", alpha = 0.4, bins = 30) +
            ggplot2::ggtitle(paste("AUCinf, model =", this_model))

          ggplot2::ggsave(file.path(run_dir, paste0("Model_", this_model, "_AUCinf_histogram_by_treatment.jpeg")), AUCinf_plot, device = "jpeg", width = 9, height = 6)
          AUClast_plot <- ggplot2::ggplot(all_NCAs_this_model, ggplot2::aes(x = AUClast), ) +
            ggplot2::geom_histogram(ggplot2::aes(color = treatment, fill = treatment), position = "identity", alpha = 0.4, bins = 30) +
            ggplot2::ggtitle(paste("AUClast, model =", this_model))

          ggplot2::ggsave(file.path(run_dir, paste0("Model_", this_model, "_AUClast_histogram_by_treatment.jpeg")), AUClast_plot, device = "jpeg", width = 9, height = 6)
        }
      }
      Cmax_plot <- ggplot2::ggplot(all_NCAs, ggplot2::aes(x = Cmax), ) +
        ggplot2::geom_histogram(ggplot2::aes(color = treatment, fill = treatment), position = "identity", alpha = 0.4, bins = 30) +
        ggplot2::ggtitle("Cmax, all Models")

      ggplot2::ggsave(file.path(run_dir, "All_models_Cmax_histogram_by_treatment.jpeg"), Cmax_plot, device = "jpeg", width = 9, height = 6)
      AUCinf_plot <- ggplot2::ggplot(all_NCAs, ggplot2::aes(x = AUCinf), ) +
        ggplot2::geom_histogram(ggplot2::aes(color = treatment, fill = treatment), position = "identity", alpha = 0.4, bins = 30) +
        ggplot2::ggtitle("AUCinf, all Models")

      ggplot2::ggsave(file.path(run_dir, "All_models_AUCinf_histogram_by_treatment.jpeg"), AUCinf_plot, device = "jpeg", width = 9, height = 6)
      AUClast_plot <- ggplot2::ggplot(all_NCAs, ggplot2::aes(x = AUClast), ) +
        ggplot2::geom_histogram(ggplot2::aes(color = treatment, fill = treatment), position = "identity", alpha = 0.4, bins = 30) +
        ggplot2::ggtitle("AUClast, all Models")

      ggplot2::ggsave(file.path(run_dir, "All_models_AUClast_histogram_by_treatment.jpeg"), AUClast_plot, device = "jpeg", width = 9, height = 6)
    },
    error = function(cond) {
      message("Failed in make_NCA_plots, error message = ", cond)
    }
  )
  return()
}

#' run_mbbe_json
#'
#' Runs MBBE from a json file of options
#' Calls run_mbbe
#'
#' @param Args.json, path to JSON file with arguments
#' @export
#' @details
#' json file structure is e.g., \cr
#' {\cr
#' "run_dir": "c:/fda/mbbe",\cr
#' "model_source": ["U:/fda/mbbe/mbbe/inst/examples/NM_05D01_11.mod",\cr
#'                  "U:/fda/mbbe/mbbe/inst/examples/NM_05D01_05.mod",\cr
#'                  "U:/fda/mbbe/mbbe/inst/examples/NM_04_085.mod",\cr
#'                  "U:/fda/mbbe/mbbe/inst/examples/NM_05D01_12.mod",\cr
#'                  "U:/fda/mbbe/mbbe/inst/examples/NM_04_090.mod"],\cr
#' "num_parallel":       32,\cr
#' "crash_value":   999999,\cr
#' "nmfe_path": "c:/nm74g64/util/nmfe74.bat",\cr
#' "delta_parms":      0.2,\cr
#' "use_check_identifiable": true,\cr
#' "NCA_end_time":       72,\cr
#' "rndseed":        1,\cr
#' "use_simulation_data": true,\cr
#' "simulation_data_path": "U:/fda/mbbe/mbbe/inst/examples/data_sim.csv",\cr
#' "ngroups":        4,\cr
#' "samp_size":      100,\cr
#' "reference_groups": [ 1, 2 ],\cr
#' "test_groups": [3, 4 ],\cr
#' "plan": "multisession",\cr
#' "alpha_error": 0.05,\cr
#' "NTID": false,\cr
#' "save_output": true,\cr
#' "model_averaging_by": "study"\cr
#' }
#'
#' @examples
#' \dontrun{
#' run_mbbe_json(Args.json)
#' }
run_mbbe_json <- function(Args.json) {
  if (!file.exists(Args.json)) {
    message(Args.json, " file not found, exiting")
  } else {
    Args <- RJSONIO::fromJSON(Args.json)
    all_args <- list(
      Args$crash_value, Args$ngroups, Args$reference_groups, Args$test_groups, Args$num_parallel, Args$samp_size, Args$run_dir, Args$model_source,
      Args$nmfe_path, Args$delta_parms, Args$use_check_identifiable, Args$NCA_end_time, Args$rndseed, Args$use_simulation_data,
      Args$plan, Args$alpha_error, Args$NTID, Args$model_averaging_by, Args$save_output
    )
    if (any(sapply(all_args, is.null))) {
      message(
        "required values in json file are: \n",
        "run_dir",
        "model_source",
        "num_parallel",
        "crash_value",
        "nmfe_path",
        "delta_parms",
        "use_check_identifiable",
        "NCA_end_time",
        "rndseed",
        "use_simulation_data",
        "ngroups",
        "samp_size",
        "reference_groups",
        "test_groups",
        "plan",
        "alpha_error",
        "NTID",
        "Args$model_averaging_by",
        "save_output and ",
        "simulation_data_path, if use_simulation_data is true\nsee documentation for details\n exiting"
      )
    } else {
      if (Args$use_simulation_data & is.null(Args$simulation_data_path)) {
        message("If use_simulation_data is true, simulation_data_path must be provided in json file")
      } else {
        return <- run_mbbe(
          Args$crash_value, Args$ngroups, Args$reference_groups, Args$test_groups, Args$num_parallel, Args$samp_size, Args$run_dir, Args$model_source,
          Args$nmfe_path, Args$delta_parms, Args$use_check_identifiable, Args$NCA_end_time, Args$rndseed, Args$use_simulation_data, Args$simulation_data_path,
          Args$plan, Args$alpha_error, Args$NTID, Args$model_averaging_by, Args$save_output
        )
        file_out <- data.frame(return$Cmax_power, return$AUClast_power, return$AUCinf_power)
        colnames(file_out) <- c("Cmax power", "AUClast power", "AUCinf power")
        write.csv(file_out, file.path(Args$run_dir, paste0("MBBEpower", stringr::str_replace_all(Sys.time(), ":", "-"), ".csv")), row.names = FALSE, quote = FALSE)
        message(paste0(capture.output(return), collapse = "\n"))
      }
    }
  }
}

#' run_mbbe
#'
#' Runs MBBE, typically called by run_mbbe_json, using a json file that includes
#' the required options
#' @param crash_value - value to be returned for BIC in models that crash, either in bootstrap or simulation
#' @param nmodels - number of models for model averaging
#' @param ngroups - number of groups in simulated data, e.g., ABBA design has 4 groups
#' @param reference_groups - which of the groups are reference formulation e.g., ABBA design might be c(2,3)
#' @param test_groups - which of the groups are test formulation e.g., ABBA design might be c(1,4)
#' @param num_parallel - number of NONMEM (bootstrap and simulation) to be run in parallel
#' @param samp_size - size of bootstrap and simulation sample
#' @param run_dir - where to run NONMEM
#' @param model_source - NONMEM control files to be used/averaged over
#' @param nmfe_path - path to nmfe??.bat
#' @param delta_parms - difference in parameters used to define identifiability
#' @param use_check_identifiable - logical, whether to use identifiability (defined by Aoki (https://www.page-meeting.org/default.asp?abstract=5951))
#' @param NCA_end_time - NCA calculation will start at 0, end at NCA_end_time
#' @param rndseed - random seed
#' @param use_simulation_data - logical, whether to used (True) a seperate data set for trial simulation or (False) the bootstrap data set
#' @param simulation_data_path - path to simulation data set
#'
#'
#' @return
#' Cmax_power,
#' AUClast_power,
#' AUCinf_power,
#' run_dir,
#' Num_identifiable,
#' BICS
#' @export
#'
run_mbbe <- function(crash_value, ngroups,
                     reference_groups, test_groups,
                     num_parallel, samp_size,
                     run_dir, model_source,
                     nmfe_path, delta_parms,
                     use_check_identifiable,
                     NCA_end_time, rndseed,
                     use_simulation_data, simulation_data_path,
                     plan = c("multisession", "sequential", "multicore"),
                     alpha_error = 0.05, NTID,
                     model_averaging_by, save_output = FALSE) {
  tictoc::tic()
  if (unlink(run_dir, recursive = TRUE, force = TRUE)) {
    stop("Unable to delete run directory ", run_dir)
  }

  # setup all models
  if (plan == "multisession") {
    old_plan <-
      future::plan(future::multisession, workers = num_parallel)
  } else if (plan == "multicore") {
    old_plan <- future::plan(future::multicore, workers = num_parallel)
  } else if (plan == "sequential") {
    old_plan <- future::plan(future::sequential)
  }

  oldOptions <- options()
  on.exit(options(oldOptions))
  options(future.globals.onReference = "error")
  on.exit(future::plan(old_plan))
  message(format(Sys.time(), digits = 0), " Start time\nModel file(s) = ", toString(model_source), "\nreference groups = ", toString(reference_groups), "\ntest groups = ", toString(test_groups))
  if (!model_averaging_by %in% c("study", "subject")) {
    stop("Error, model_averaging is ",model_averaging_by, " model_averaging_by must be one of study or subject, exiting")
  } else {
    message("Model averaging will be by ", model_averaging_by)
  }
  message("Bootstrap/Monte Carlo sample size = ", samp_size, "\nnmfe??.bat path = ", nmfe_path, "\nUse_check_identifiability = ", use_check_identifiable)
  message("Narrow Therapeutic Index  = ", NTID)
  message("Alpha error rate for bioqulvalence testing = ", alpha_error)
  message("Number parallel runs for bootstrap, simulations and NCA ", num_parallel)
  if (use_check_identifiable) {
    message("Delta parameter for use_check_identifiable = ", delta_parms)
  }
  if (use_simulation_data) {
    message("Simulation data set = ", simulation_data_path)
  } else {
    message("Original analysis data set will be used for simulation")
  }

  path_parents <- split_path(run_dir)
  cur_path <- path_parents[1]
  for (this_parent in 2:length(path_parents)) {
    cur_path <- file.path(cur_path, path_parents[this_parent])
    if (!dir.exists(cur_path)) {
      dir.create(cur_path)
    }
  }

  set.seed(rndseed)

  msg <- check_requirements(
    run_dir, samp_size, model_source, ngroups,
    reference_groups, test_groups, nmfe_path,
    use_check_identifiable, use_simulation_data,
    simulation_data_path
  )
  if (msg$rval) {
    message(
      "Passed requirements check\nCopying source control files from ", toString(model_source), " to ", file.path(run_dir, "modelN"),
      "\n where N is the model number"
    )
    nmodels <- copy_model_files(model_source, run_dir)

    if (nmodels > 0) {
      passes_identifiable <- c(rep(0, nmodels))
      names(passes_identifiable) <- paste0("Model", seq(1:nmodels))
      message(format(Sys.time(), digits = 0), " Sampling data 1-", samp_size, " writing data to ", file.path(run_dir, "data_sampM.csv"), " where M is the bootstrap sample number")
      sample_data(run_dir, nmodels, samp_size)
      message(format(Sys.time(), digits = 0), " Starting bootstrap runs 1-", samp_size, " in ", file.path(run_dir, "modelN", "M"), " where N is the model number and M is the sample number")
      if (!run_any_models(nmfe_path, run_dir, nmodels, samp_size, TRUE)) {
        message("Failed bootstrap")
      } else {
        # need to wait until all are done, this returns when all are started.

        message("\n", format(Sys.time(), digits = 0), " Getting bootstrap model parameters, samples 1-", samp_size)
        parms <- get_parameters(
          run_dir, nmodels,
          samp_size, delta_parms,
          crash_value, use_check_identifiable,
          passes_identifiable
        )
        base_models <- get_base_model(run_dir, nmodels) # get all nmodels base model
        message(format(Sys.time(), digits = 0), " Constructing simulation  models in ", file.path(run_dir, "MBBEsimM"), " where M is the simulation number")
        final_models <- write_sim_controls(run_dir, parms, base_models, samp_size, use_simulation_data, simulation_data_path) # don't really do anything with final models, already written to disc
        message(format(Sys.time(), digits = 0), " Running simulation models 1-", samp_size, " in ", file.path(run_dir, "MBBEsimM"), " where M is the simulation number")
        if (run_any_models(nmfe_path, run_dir, NULL, samp_size, FALSE, plan, num_parallel)) {
          message("\n", format(Sys.time(), digits = 0), " Calculating NCA parameters for simulations 1-", samp_size, ", writing to ", file.path(run_dir, "MBBEsimM", "NCAresultsM.csv"), ",  where M is the simulation number")
          if (model_averaging_by == "subject") {
            message("\n Note the NCA parameters are by study, prior random selection by subject for model averaging")
          }

          # note would like to do NCA parallel, so don't end
          calc_NCA(
            run_dir,
            ngroups,
            reference_groups,
            test_groups,
            NCA_end_time,
            samp_size, plan,
            num_parallel
          )
          message(format(Sys.time(), digits = 0), " Done calculating NCA parameters, making plots")

          make_NCA_plots(
            parms$BICS,
            run_dir,
            samp_size,
            nmodels,
            reference_groups,
            test_groups
          )

          message(format(Sys.time(), digits = 0), " Done making plots, calculation power")
          # read NCA output and do stats
          all_results <- calc_power(run_dir, samp_size,
            alpha = alpha_error,
            model_averaging_by, NTID = NTID
          )
        } else {
          all_results <- NULL
        }
        future::plan(old_plan)
        if (!is.null(all_results)) {
          output_file <- file.path(run_dir, "All_results.csv")
          count <- 0
          while (file.exists(output_file) & count < 20) {
            count <- count + 1
            file.remove(output_file)
            Sys.sleep(0.25)
          }
          count <- 0
          while (!file.exists(output_file) & count < 20) {
            write.csv(all_results, file = output_file, quote = FALSE)
            count <- count + 1
            Sys.sleep(0.25)
          }
          if (!file.exists(output_file)) {
            message("Unable to write to ", output_file)
          }
          Cmax_power <- all_results %>%
            dplyr::filter(Cmax_BE != -999) %>%
            dplyr::summarise(Power = mean(Cmax_BE))
          AUClast_power <- all_results %>%
            dplyr::filter(AUClast_BE != -999) %>%
            dplyr::summarise(Power = mean(AUClast_BE))
          AUCinf_power <- all_results %>%
            dplyr::filter(AUCinf_BE != -999) %>%
            dplyr::summarise(Power = mean(AUCinf_BE))
          power <- c(Cmax_power, AUCinf_power, AUClast_power)
          write.csv(power, file.path(run_dir, "Power.csv"))
        } else {
          message("Failed in calc_power")
          message(msg$msg, " exiting")
          Cmax_power <- NULL
          AUClast_power <- NULL
          AUCinf_power <- NULL
        }
      }
    }
  } else {
    message(msg$msg, " exiting")
    Cmax_power <- NULL
    AUClast_power <- NULL
    AUCinf_power <- NULL
  }
  tictoc::toc()
  message("Done at ", Sys.time())
  return(
    list(
      "Cmax_power" = Cmax_power,
      "AUClast_power" = AUClast_power,
      "AUCinf_power" = AUCinf_power,
      "run_dir" = run_dir,
      "Num_identifiable" = parms$passes_identifiable,
      "BICS" = parms$BICS
    )
  )

}
