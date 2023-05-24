# The term 'model' refers the the source model, the models to be averaged. The term sample refers to the bootstrap sample of the data set
# and the model run against that bootstrap sample of the data set, and the set of Monte Carlo simulations based on the best model from the
rm(list = ls())
library(RJSONIO)  # fromJSON
library(dplyr)  # all
library(stringr)  # str_trim, str_replace
library(xml2)
library(future)  # future
library(PKNCA)
library(DescTools)  # splitPath
library(installr)  # is.empty

delete_files <-  function(folder){
  # need .ext for use_identifiable, keep .lst, .xml, .mod, GMSF
  # but can't delete all, need to keep $TABLE
  try({
  delfiles <- c("FDATA","FCON","FSUBS","FSUBS.o","FSIZES", "FREPORT","FSTREAM", "GFCOMPILE.BAT", "INTER","nmprd4p.mod" ) 
  delfiles <- c(delfiles, dir(path = folder, pattern = "*.exe"), dir(path = folder, pattern = "*.f90"), dir(path = folder, pattern = "*.grd"), 
                dir(path = folder, pattern="*.shk"), dir(path = folder, pattern="*.cpu"), dir(path = folder, pattern="*.shm"), 
                dir(path = folder, pattern="*.lnk"), dir(path = folder, pattern = "*.phi"))
  file.remove(file.path(folder, delfiles))
  if(dir.exists(file.path(folder,"temp_dir"))){
    unlink(file.path(folder,"temp_dir"),recursive = TRUE, force = TRUE)
  }
  })
}
#' get_block
#' Get a block (e.g., $DATA) from a control file with multiple lines
#' return block as a single character string
#'
#' @param stem Block title, e.g., $DATA to look for
#' @param control control file text 
#'s
#' @examples
#' get_block('$DATA',c('$PROB test','$INPUT ...)) 
get_block <- function(stem, control) {

    tryCatch({
        nlines <- length(control)
        for (this_line in 1:nlines) {
            # remove comments
            comment.pos <- gregexpr(";", control[this_line])[[1]][1]
            if (comment.pos > 0) {
                control[this_line] <- stringr::str_trim(substr(control[this_line], 1, comment.pos - 1))
            } else {
                control[this_line] <- stringr::str_trim(control[this_line])
            }
        }
        start_line <- grep(paste0("^\\", stem), control)
        next_line <- grep("^\\$", control[start_line + 1:length(control)])[1]
        if (installr::is.empty(next_line)) {
            # $stem is last
            last_line <- length(control)
        } else {
            last_line <- start_line + next_line
        }
        block <- ""
        for (line in start_line:(last_line - 1)) {
            block <- paste(block, control[line])
        }
        return(block)
    }, error = function(err) {
        return(FALSE)
    })
}
#' check_requirements
#' 1. Is there exactly 1 .mod file in each source folder
#' 2. Does that .mod file contain ;;;;.*Start EST
#' 3. is the sum of reference_groups and test_groups == ngroups?
#' 4. Any duplicates in Reference and Test groups?
#' 5. check if data file is present
#' 6.  check if nmfe path is correct
#' 7. check for saddle_reset if requested
#' 8. check for repeat IDs in data set 
#' 9. if use_simulation_data, see if data available
#' 
#' @param source.dir source  model files
#' @param nmodels number of models
#' @param ngroups number of groups in data set
#' @param reference_groups Which groups are reference
#' @param test_groups Which groups are Test
#' @param nmfe_path nmfe??.bat
#' @param use_check_identifiable - logical, whether to check for identifability, will check if SADDLE_RESET is in $EST
#' @param use_simulation_data logical, if the simulation will be done with a different data set than the bootstrap
#' @param simulation_data_path, if use_simulation_data, this is the path to the data file
#' @examples
#' check_requirements('c:/models',5,'c:/nm744/util',TRUE)) 
check_requirements <- function(model_list, ngroups, reference_groups, test_groups, nmfe_path, use_check_identifiable, use_simulation_data, simulation_data_path = NULL) {

    msg <- list()
    result = tryCatch({
        if(is.null(model_list)){
          return(list(rval = FALSE, msg = paste("Model list is NULL, error in json file?, exiting")))
        }
        if (!file.exists(nmfe_path)) {
            # if DOS path, convert to R/linux
            return(list(rval = FALSE, msg = paste("Cannot find nmfe?? at", nmfe_path, ", exiting")))
        }
        # check number in Reference and Test groups
        if (sum(length(reference_groups), length(test_groups)) != ngroups) {

            return(list(rval = FALSE, msg = paste("number of Reference groups ", length(reference_groups), "+ Test groups ", length(test_groups),
                "doesn't equal the number of groups", ngroups, ", exiting")))
        }
        # no duplicated in Reference and Test groups
        if (anyDuplicated(c(reference_groups, test_groups)) > 0) {
            return(list(rval = FALSE, msg = paste("There are duplicated group numbers between Reference and Test group, exiting")))
        }
        if (anyDuplicated(reference_groups) > 0) {
            return(list(rval = FALSE, msg = paste("There are duplicated group numbers in the Reference group, exiting")))
        }
        if (anyDuplicated(test_groups) > 0) {
            return(list(rval = FALSE, msg = paste("There are duplicated group numbers in the test group, exiting")))
        }
    }, error = function(err) {

        msg = list(rval = FALSE, msg = paste("Error in finding nfme??.bat,", nmfe_path, ", exiting"))
    })
    nmodels <- length(model_list)
    for (this_model in 1:nmodels) {
        # need to fix this for modellist
        result <- tryCatch({
            for (this_model in 1:nmodels) {
                if (!file.exists(model_list[this_model])) {
                  return(list(rval = FALSE, msg = paste("Cannot find", file.path(source.dir, paste0("model", this_model)), ", exiting")))
                } else {

                  con <- file(model_list[this_model], "r")
                  suppressWarnings(control <- readLines(con, encoding = "UTF-8"))
                  close(con)
                  data_line <- get_block("$DATA", control)
                  data_line <- stringr::str_trim(stringr::str_replace(data_line, "\\$DATA", ""), side = "both")
                  any.quotes = grep("^\"", data_line)
                  if (length(any.quotes) > 0) {
                    # find 2nd
                    pos <- gregexpr(pattern = "\"", data_line)
                    data_file <- substr(data_line, 1, pos[[1]][2])
                  } else {
                    # find first white space
                    pos <- gregexpr(pattern = "\\s", data_line)
                    data_file <- stringr::str_trim(substr(data_line, 1, pos[[1]][1]), side = "both")
                  }
                  if (!file.exists(data_file)) {
                    msg <- append(msg, list(rval = FALSE, msg = paste0("Cannot find ", data_file, ", exiting")))
                  } else {
                    # repeat IDs??

                    data <- read.csv(data_file, stringsAsFactors = FALSE)
                    newIDs <- data %>%
                      dplyr::select(ID) %>%
                      mutate(newID = if_else(ID == lag(ID), 0, 1)) %>%
                      mutate(newID = if_else(is.na(newID), 1, newID)) %>%
                      dplyr::filter(newID == 1)
                    num_newIDs <- dim(newIDs)[1]
                    num_oldIDs <- dim(data %>%
                      dplyr::distinct(ID))[1]
                    if (num_newIDs != num_oldIDs) {
                      msg = append(msg, list(rval = FALSE, msg = paste0("There appears to be repeat IDs in", data_file, "\n, for bootstrap sampling IDs must not repeat, exiting")))
                    }
                  }

                  # $EST must be on one line, need to fix this contains ;;;; Start EST?
                  contains_start <- grepl(";;;;.*Start\\s+EST", control, ignore.case = TRUE)
                  if (!any(contains_start)) {
                    msg = append(msg, list(rval = FALSE, msg = paste0("The control file ", model_list[this_model], " does not contain \";;;; Start EST\", required before the $EST and the $THETA records and after $OMEGA and $SIGMA")))
                  }
                  if (use_check_identifiable) {
                    EST.line <- get_block("$EST", control)
                    if (!grepl("SADDLE_RESET\\s*=\\s*1", EST.line, fixed = FALSE)) {
                      msg = append(msg, list(rval = FALSE, msg = paste0("Identifiability check requested, but SADDLE_RESET not set to 1 in ",
                        model_list[this_model], ", exiting")))
                    }
                  }
                  if (use_simulation_data) {
                    if (!file.exists(simulation_data_path)) {
                      msg <- append(msg, list(rval = FALSE, msg = paste0("Cannot find simulation data file ", simulation_data_path, ", exiting")))
                    }
                  }
                }
            }
        }, error = function(err) {
            msg = append(msg, list(rval = FALSE, msg = paste("Error in number of model files in", model_list, ", exiting")))
        })
        if (length(msg) > 0) {
            return(msg)
        } else {
            return(list(rval = TRUE, msg = paste("passes requirements check")))
        }
    }
}
#' split_path
#' split a path into parents
#' 
#' @param  path, string, path name to be split
#' @param mustWork logical
#' @returns list of path parents
#' @examples
split_path <- function(path, mustWork = FALSE) {
    output <- c(strsplit(dirname(normalizePath(path, mustWork = FALSE)), "/|\\\\")[[1]], basename(path))
    return(output)
}
#' copy_model_files
#' copy NONMEM model files from a source to a run directory
#' Need to look at this for OS portability??
#' @param  model_source folder name were subfolders (model1..modelnmodels) with model files can be found. There must be exactly one .mod file (with any stem) in the folder
#' @param run_dir Folder where models are to be run
#' @return nmodels how many models are there
#' @examples
#' copy_model_files('c:/modelaveraging','c:/modelaveraging/run',4)
copy_model_files <- function(model_source, run_dir) {

    if (!file.exists(run_dir)) {
        # split run_dir, check each parent
        normpath <- SplitPath(run_dir)$normpath
        parents <- split_path(normpath)
        nparents <- length(parents)
        cur_path <- parents[1]  # should be the drive
        for (this.parent in 2:nparents) {
            cur_path = file.path(cur_path, parents[this.parent])
            if (!file.exists(cur_path)) {
                dir.create(cur_path)
            }
        }
    }
    nmodels <- length(model_source)
    for (this_model in 1:nmodels) {
        bs_dir <- file.path(run_dir, paste0("model", this_model))  # bootstrap directory
        if (dir.exists(bs_dir)) {
            unlink(bs_dir, recursive = TRUE, force = TRUE)
        }
        dir.create(bs_dir)
        modfile <- model_source[this_model]
        file.copy(modfile, file.path(bs_dir, paste0("bs", this_model, ".mod")))
    }
    return(nmodels)
}

# copy_model_files <- function(model_source,run_dir,nmodels){ if(!file.exists(run_dir)){ # split run_dir, check each parent normpath =
# SplitPath(run_dir)$normpath parents = split_path(normpath) nparents = length(parents) cur_path = parents[1] for(this.parent in
# 2:nparents){ cur_path = file.path(cur_path,parents[this.parent]) if(!file.exists(cur_path)){ dir.create(cur_path) } } } for(this_model
# in 1:nmodels){ sim_dir <- file.path(run_dir,paste0('model',this_model)) if(dir.exists(sim_dir)){ unlink(sim_dir,recursive = TRUE, force
# = TRUE) } dir.create(sim_dir) source.dir <- file.path(model_source,paste0('model',this_model)) modfile <-
# file.path(model_source,paste0('model',this_model),list.files(source.dir))
# file.copy(modfile,file.path(sim_dir,paste0('bs',this_model,'.mod'))) } }

#' create bootstrap samples of NONMEM data set, placed in run_dir, file names = data_sampN.csv
#' 
#' @param run_dir Folder where models are to be run
#' @param samp_size how many bootstrap samples, default is 100
#' @param nmodels how many models are there
#' @examples
#' sample_data('c:/modelaveraging',100,4)
sample_data <- function(run_dir, nmodels, samp_size) {

    con <- file(file.path(run_dir, "model1", "bs1.mod"), "r")
    suppressWarnings(control <- readLines(con, encoding = "UTF-8"))
    close(con)
    data_line <- stringr::str_trim(control[grep("\\$DATA", control)])
    data_line <- stringr::str_trim(stringr::str_replace(data_line, "\\$DATA", ""), side = "both")
    # if file name is quoted, just put out part in in quotes, otherwise get first white space
    any.quotes = grep("^\"", data_line)
    if (length(any.quotes) > 0) {
        # find 2nd
        pos <- gregexpr(pattern = "\"", data_line)
        data_file <- substr(data_line, 1, pos[[1]][2])
        rest.of.data_line <- substr(data_line, pos[[1]][2], nchar(data_line))
    } else {
        # find first white space
        pos <- gregexpr(pattern = "\\s", data_line)
        data_file <- stringr::str_trim(substr(data_line, 1, pos[[1]][1]), side = "both")
        rest.of.data_line = substr(data_line, pos[[1]][1], nchar(data_line))
    }
    # get datafile possibly different data files in different models? not supported at this time
    org.data <- read.csv(data_file, header = TRUE, stringsAsFactors = FALSE)
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
            con <- file(file.path(run_dir, paste0("model", this_model), paste0("bs", this_model, ".mod")), "r")
            suppressWarnings(control <- readLines(con, encoding = "UTF-8"))
            close(con)
            data_line <- grep("\\$DATA", control)
            # get rest of data line !!! note, will need to read entire $DATA block, maybe on more than one line!!!!
            newdata_line <- paste0("$DATA ", file.path("..", "..", paste0("data_samp", this_samp, ".csv ", rest.of.data_line)))
            control[data_line] <- newdata_line
            # remove covariance
            cov_line <- grep("\\$COV", control)
            if (length(cov_line) > 0)
                control[cov_line] = ";; no covariance step"
            # and any tables
            control <- stringr::str_replace(control, "$TABLE", ";$TABLE")
            dir.create(file.path(run_dir, paste0("model", this_model), this_samp))

            con <- file(file.path(run_dir, paste0("model", this_model), this_samp, paste0("bsSamp", this_samp, ".mod")), "w")
            writeLines(control, con)
            close(con)
        }
    }
}

#' run the bootstrap models/samples
#' 
#' @param nmfe_path path to nmfe??.bat
#' @param run_dir Folder where models are to be run 
#' @param nmodels how many models are there
#' @param samp_size how many samples are there
#' @param numParallel how many to run in parallel, default = availableCores()
#' @examples
#' run.bootstrap('c:/nmfe744/util/nmfe74.bat', 'c:/modelaveraging',8)
run.bootstrap <- function(nmfe_path, run_dir, nmodels, samp_size, numParallel = availableCores()) {
  run_one_model <- function(run_dir, nmfe_path,this_model, this_samp){
    setwd(file.path(run_dir, paste0("model", this_model), this_samp))
    command <- paste0(nmfe_path, " bsSamp", this_samp, ".mod bsSamp", this_samp, ".lst -nmexec=", paste0("bsSamp", this_model,
                                                                                                         "_", this_samp, ".exe"))
    shell(command, wait = TRUE)
    delete_files(getwd())
  }
  
    run_count <- 0
    rval <- tryCatch({
        plan(multisession, workers = numParallel)
        for (this_model in 1:nmodels) {
            pb <- txtProgressBar(min = 0, max = samp_size, initial = 0)
            message(Sys.time()," Starting first bootstrap sample model ", this_model)
            for (this_samp in 1:samp_size) {
              run_count <- run_count + 1
                future::future({
                  futs[run_count] <- run_one_model (run_dir, nmfe_path, this_model, this_samp)
                })
                setTxtProgressBar(pb, this_samp)
            }
            close(pb)
            message(Sys.time(), " Starting last bootstrap sample model ", this_model) 
        } 
        close(pb)
        return(TRUE)
    }, error = function(cond) {
        return(FALSE)
    })
    
    # wait until all done
    for(this_fut in futs){
      while (!resolved(this_fut)) { 
        Sys.sleep(1) 
      }
    }
    
    plan(sequential)
}


#' get_parameters
#' read the xml file for each bootstrap/model, calculate BIC
#'  return list of BICs (with the BEST column) and parameters for all models/samples
#'  
#' @param run_dir Folder where models are to be run 
#' @param nmodels how many models are there
#' @param samp_size how many samples are there 
#' @param delta_parms criteria for identifability test
#' @param crash_value value for failed calculation
#' @examples
#' get_parameters('c:/modelaveraging',4,100,0.1)
get_parameters <- function(run_dir, nmodels, samp_size, delta_parms, crash_value, use_check_identifiable) {

    BICS <- data.frame(matrix(crash_value, nrow = samp_size, ncol = nmodels + 3))
    colnames(BICS) <- c(paste0("Model", seq(1:nmodels)), "Best","Max_Delta_parm","Max_Delta")
    # BIC=k*ln(n) -2LL where k is is the number of estimated parameters and n is the number of data
    nparms <- data.frame(ntheta = as.integer())

    # only need parameters for best model, but need to read xml file to get the number of observtions and parameter, so may as well get
    # THETA etc now for all
    nfailed_ident <- 0  # number of samples in this model that fail the identifabilit test
    # do by sample first only save parameters for the best model
    selected_parameters <- list()  # parameters for selected model, 
    this_samp = 1
    for (this_samp in 1:samp_size) {
        num_successful <- 1  # only save parameters if model finished, doesn't need to converge, just finish, this is the number of successful samples for this model
        parameters_this_sample <- list()
        this_model = 1
        all_identifiables <- rep(TRUE, nmodels)
        for (this_model in 1:nmodels) {
           
            xml_file <- file.path(run_dir, paste0("model", this_model), this_samp, paste0("bsSamp", this_samp, ".xml"))
            count <- 0
            while (!file.exists(xml_file) & count < 200) {
               Sys.sleep(1)
               count <- count + 1 
             }
            if(!file.exists(xml_file)){
              BICS[this_samp, this_model] <- crash_value
              next
            }
            rval <- tryCatch({
              count <- 0 
              boolFalse <- FALSE
              # wait for files to close?? otherwise get Error in file(file, "rt") : cannot open the connection
              while(boolFalse == FALSE & count < 50) {
                data <-  NULL  
                count <- count + 1
                tryCatch({
                  data <- xml2::read_xml(xml_file, encoding = "ASCII")
                  boolFalse <- TRUE
                }, error = function(e){
                     Sys.sleep(1) 
                }, finally = {}) # fails will throw error below in outer tryCatch
              }
               
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
                  nparms <- rbind(nparms, data.frame(ntheta = ntheta)) 
              }
              parameters_this_model <- data.frame(matrix(-9999, ncol = ntheta))
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
              }else{ 
                # leave as -9999 from above parameters_this_sample[[this_model]] <- -9999
              }
  
  
              # if using saddle_reset, test for identifiability
              if (!use_check_identifiable) {
                  identifiable_ok <- TRUE
              } else {
                  # run_dir,this_model,this_sample,delta_parms
                  identifiable_ok <- check_identifiable(run_dir, this_model, this_samp, delta_parms, ntheta)
                  all_identifiables[this_model] <- identifiable_ok["passes"]
                  if (!identifiable_ok$passes) {
                    # first just get the number of parameters and observations for this model
                    BICS[this_samp, this_model] <- crash_value
                  }
            }
              },error=function(e){
                all_identifiables[this_model] <- identifiable_ok["passes"]
                BICS[this_samp, this_model] <- crash_value
              
            })
             
        }
        # and select best model, -9999 if all crashed
        # seems to have copied model 8 parms from model 7
        if (all(BICS[this_samp, 1:nmodels] == crash_value)) {
            BICS$Best[this_samp] <- -9999 
            warning("No models successful for model ", this_samp, " identifiability = ", toString(all_identifiables))
            selected_parameters[[this_samp]] <- -9999
        } else {
            best <- which.min(BICS[this_samp, 1:nmodels])
            BICS$Best[this_samp] <- best
            selected_parameters[[this_samp]] <- unlist(parameters_this_sample[best])
        }
        BICS$Max_Delta[this_samp] <- identifiable_ok$max_delta
        BICS$Max_Delta_parm[this_samp] <- identifiable_ok$max_delta_parm

    }
    # write out results, 
    write.csv(BICS, file.path(run_dir, "BICS.csv"), quote = FALSE)
    
    # write out parameters, can't do write.csv, as rows have different number of parameters
    conn <- file(file.path(run_dir, "Parameters.csv"), open="w")
    
    newlist <- lapply(seq_len(length(selected_parameters)), function(i){
      
      temp = lapply(seq_len(length(selected_parameters[[i]])), function(j) {
        temp <- c(selected_parameters[[i]][[j]])
      })
      
      writeLines(text=paste(i,paste(temp, collapse=","),sep=","), con=conn, sep="\n")  
    })
    close(conn)
    rval <- list(BICS = BICS, parameters = selected_parameters)
    return(rval)
}

#' for each model get the control file used for the bootstrap
#' return a list of the models  
#' @param nmodels how many models are there 
#' @examples
#' get_base_model(4)
get_base_model <- function(nmodels) {

    # need error trapping for no ;;;;;.*Start EST
    base_models <- vector(mode = "list", length = 0)  # all but $THETA, $SIM, $TABLE 

    for (this_model in 1:nmodels) {
      tryCatch({
        con <- file(file.path(paste0("model", this_model), paste0("bs", this_model, ".mod")), "r")
        suppressWarnings(control <- readLines(con, encoding = "UTF-8"))
        control <- control[1:grep(";;;;.*Start EST", control)]
        close(con)
        base_models <- append(base_models, list(control))
      },error = function(e){
        stop("Cannot find ", file.path(paste0("model", this_model), paste0("bs", this_model, ".mod")), "exiting")
      })
    }
    return(base_models)


}

#' Edits the best based model for each sample, replaces the original parameters with the bootstrap parameters
#' and adds the $SIM and $TABLE
#'  
#' @param run_dir Folder where models are to be run 
#' @param parms list that include BICs and parameters for each sample/model
#' @param base_models list of the text in each model used for the bootstrap
#' @param samp_size how many samples are there 
#' @examples
#' write_sim_controls('c:/modelaveraging',parms,base_models,100)
write_sim_controls <- function(run_dir, parms, base_models, samp_size, use_simulation_data, simulation_data_path = NULL) {

    nmodels <- length(base_models)
    final_models <- vector(mode = "list", length = samp_size)
    model.indices <- rep(0, nmodels)  # which parameter set to use when this model is selected, roll over if not enough samples,
    # shouldn't happen, as we shouldn't use all parameter sets for any model

    for (this_samp in 1:samp_size) {
        which.model <- parms$BICS$Best[this_samp]
        if (which.model == -9999) {
            message("skipping Sample #", this_samp, " no acceptable model found")
        } else {
            # use BS parameters, when you run out (as some BS samples fail), just start over, so need different random seed in $SIM model
            # is first index, e.g., parameters[[2]]$THETA1[1] is THETA[1] for model 2, first sample
            model.indices[which.model] <- model.indices[which.model] + 1
            # if not enough model parameters (because some model crashed?), recycle. But, shouldn't need to as model that crashes should
            # have crash_value BIC
            if (model.indices[which.model] > length(parms$parameters[[1]])[1])
                model.indices[which.model] <- 1
            full_control <- base_models[which.model][[1]]
            # need to get sequence in here, calculated in $PK
            seed <- round(runif(1, 0, 10000), 0)
            full_control <- c(full_control, paste("$SIM ONLYSIM (", seed, ")"), "$THETA")
            ntheta <- length(parms$parameters[[which.model]])

            if (use_simulation_data) {
                # get $DATA
                start_line <- grep("^\\$DATA", full_control)
                next_line <- grep("^\\$", full_control[start_line + 1:length(full_control)])[1]
                if (installr::is.empty(next_line)) {
                  # $DATA is last
                  last_line <- length(full_control[[1]])
                } else {
                  last_line <- start_line + next_line
                }
                # replace with simulation data set
                full_control <- c(full_control[1:(start_line - 1)], paste("$DATA", simulation_data_path, "IGNORE=@"), full_control[(last_line):length(full_control)])
            }
            for (this_parm in 1:ntheta) {
                full_control <- c(full_control, paste0(parms$parameters[[which.model]][this_parm], "  ;; THETA(", this_parm, ")"))
            }
            full_control <- c(full_control, "$TABLE ID TIME GROUP OCC SEQ DV EVID NOPRINT NOAPPEND FILE=OUT.DAT ONEHEADER")
            sim_dir <- file.path(run_dir, paste0("Sim", this_samp))
            if (dir.exists(sim_dir)) {
                unlink(sim_dir, recursive = TRUE, force = TRUE)
            }
            if (dir.exists(sim_dir)) {
                unlink(sim_dir, recursive = TRUE, force = TRUE)
            }
            dir.create(sim_dir)
            con <- file(file.path(sim_dir, paste0("sim", this_samp, ".mod")), "w")
            writeLines(unlist(full_control), con)
            close(con)
            final_models <- append(final_models, list(full_control))
        }
    }
    return(final_models)
}

#' wait_for_bs
#' make array of all PIDs, check when all return Null 
#' Stem for bootstrap is bsSamp,  
#' Note this assume at least one run has started, could be a problem if all runs
#' are still running nmtran/compiler
#' Note that there is no time out on runs, as there models should have all run to
#' completion with the full data set, presumably they will be the bootstrap data set
#' @param nmodels number of base models
#' @param samp_size number of bootstrap sample for each model
#' @examples
#' wait.for_bs(4,100) 
# names of runs will be stem(modelnum)_(samp_num).exe
wait_for_bs = function(nmodels, samp_size) {
    total.models <- nmodels * samp_size
    models.done = 0
    pb <- txtProgressBar(min = 0, max = total.models, initial = 0)
    PID <- vector("numeric", length = nmodels * samp_size)
    this_run <- 0
    for (this_model in 1:nmodels) {
        for (this_samp in 1:samp_size) {
            this_run <- this_run + 1
            PID[this_run] <- max(0, get_pid(paste0("bsSamp", this_model, "_", this_samp, ".exe")))
        }
    }
    # loop until all are 0
    while (any(PID > 0)) {
        Sys.sleep(5)
        this_run <- 0
        for (this_model in 1:nmodels) {
            for (this_samp in 1:samp_size) {
                this_run <- this_run + 1
                if (PID[this_run] > 0) {
                  PID[this_run] <- max(0, get_pid(paste0("bsSamp", this_model, "_", this_samp, ".exe")))
                  if (PID[this_run] == 0) {
                    models.done <- models.done + 1
                    setTxtProgressBar(pb, models.done)
                  }
                }
            }
        }
    }
    close(pb)
}

#' wait_for_sim
#' Pretty much the same as wait_for_bs, except no nmodels.
#' make array of all PIDs, check when all return Null 
#' Stem for simulation is sim
#' Note this assume at least one run has started, could be a problem if all runs
#' are still running nmtran/compiler 
#' @param samp_size number of bootstrap sample for each model
#' @examples
#' wait_for_sim(100) 
# names of runs will be sim(samp_num).exe
wait_for_sim = function(samp_size) {

    PID <- vector("numeric", length = samp_size)
    this_run <- 0
    for (this_samp in 1:samp_size) {
        this_run <- this_run + 1
        PID[this_run] <- max(0, get_pid(paste0("sim", this_samp, ".exe")))
    }
    # loop until all are 0
    while (any(PID > 0)) {
        Sys.sleep(5)
        this_run <- 0
        for (this_samp in 1:samp_size) {
            this_run <- this_run + 1
            if (PID[this_run] > 0) {
                PID[this_run] <- max(0, get_pid(paste0("sim", this_samp, ".exe")))  # get_pid returns NULL is no such process

            }
        }
    }
}
#' run all simulations
#' @param nmfe_path path to nmfe??.bat
#' @param run_dir home directory
#' @param numParallel number of models to run in parallel, default = avaiableCores()
#' @param samp_size number of bootstrap sample for each model
#' @examples
#' run_simulations('c:/nmfe75/util/nmfe75.bat','c:/mbbe',8,100)
run_simulations <- function(nmfe_path, run_dir, samp_size, numParallel = availableCores()) {
  
  
  run_one_sim <- function(run_dir, nmfe_path, this_samp){
    setwd(file.path(run_dir, paste0("sim", this_samp)))
    if (file.exists(paste0("sim", this_samp, ".mod"))) {
      # no control file written if no models were acceptable
      command <- paste0(nmfe_path, " sim", this_samp, ".mod sim", this_samp, ".lst -nmexec=sim", this_samp, ".exe")
      shell(command, wait = TRUE)
      delete_files(getwd())
    }
  }
  futs <- list()
    plan(multisession, workers = numParallel)
    for (this_samp in 1:samp_size) {
        future::future({
          futs[this_samp] <- run_one_sim(run_dir, nmfe_path, this_samp)
        })
    }
    # but this returns once all nmtran runs are done, doesn't wait for nonmem run to be done
    for(this_fut in futs){
      while (!resolved(this_fut)) { 
        Sys.sleep(0.2) 
      }
    }
    plan(sequential) 
}

#' calc_NCA
#' call getNCA for each sample in 1:samp_size
#' and do NCA (Cmax, AUCin,AUClst, from 0 to end.time)
#' if check.identifability, do that with delta_parms as criteria
#' @param run_dir, string home directory  
#' @param nGroups, integer how many groups, e.g., 4 for ABBA
#' @param reference_groups, arrays which groups are reference  
#' @param test_groups, array, which groups are  test
#' @param NCA_end_time, numeric, end time for AUClast and AUCinf
#' @param samp_size,numeric integer
#' @param delta_parms, numeric absolute difference in parameters criteria for failig identifability
#' @examples
#' calc_NCA('c:/runmodels',4,c(1,2),c(3,4)),72,100,0.1)
calc_NCA <- function(run_dir, ngroups, reference_groups, test_groups, NCA_end_time, samp_size, numParallel = availableCores()) {
    pb <- txtProgressBar(min = 0, max = samp_size, initial = 0)  
  
    futs <- list() # list of futures
    plan(multisession, workers = numParallel)
    for (this_samp in 1:samp_size) {
        future::future({
            # writes to file, future can't return a value?
          futs[this_samp] <- getNCA(run_dir, this_samp, ngroups, reference_groups, test_groups, NCA_end_time)
        }) 
    }
  # wait until all done
  for(this_fut in futs){
    setTxtProgressBar(pb, this_fut)
    while (!resolved(this_fut)) { 
         Sys.sleep(0.2) 
     }
    }
    plan(sequential) 
    close(pb)
}
#' read .lst file, find out where saddle reset occurs then get parameter before reset from .ext file (Pre.parms)
#' compare with final parameters (Post.parms), see if any differ by delta_parms
#' open .lst 
#' @param run_dir home directory
#' @param this_model integer
#' @param this_sample integer
#' @param delta_parms absolute difference in parameters criteria for failing identifability
#' @examples
#' check_identifiable('c:/runmodels',1,1,0.1)
check_identifiable <- function(run_dir, this_model, this_sample, delta_parms, nparms) {
    lstfile <- file.path(run_dir, paste0("model", this_model), this_sample, paste0("bsSamp", this_sample, ".lst"))
    extfile <- file.path(run_dir, paste0("model", this_model), this_sample, paste0("bsSamp", this_sample, ".ext"))
    max_deltap <- -999
    max_delta_parmp <- -999
    if (!file.exists(lstfile) | !file.exists(extfile)) {
      message("Cannot find ", listfile, " or ", extfile, " for determining identifiability")
      return(c(passes = FALSE, max_delta = -999, max_delta_parm = -999))
    } else {
        tryCatch({
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
            ext <- read.table (extfile, header = TRUE,  skip = 1)
            Pre.parms <- ext %>%
                dplyr::filter(ITERATION == reset_iteration)
            Post.parms <- ext %>%
                dplyr::filter(ITERATION == last.iteration)
            #nparms <- length(Pre.parms)  # includes the first column - 'ITERATION'
            Passes_identifiability <- TRUE
             
            for (this_parm in 2:nparms) {
                if (Post.parms[this_parm] != 0) {
                  difference <- abs((Pre.parms[this_parm] - Post.parms[this_parm])/Post.parms[this_parm])
                  if (difference > delta_parms) {
                    Passes_identifiability <- FALSE
                  }
                  if (difference > max_deltap) {
                    max_deltap <- difference
                    max_delta_parmp <- this_parm - 1
                  }
                }
            }
            rval = c(passes = Passes_identifiability, max_delta = max_deltap, max_delta_parm = max_delta_parmp)
            names(rval ) = c("passes","max_delta","max_delta_parm")
            return(rval)
        }, error = function(cond) {
            rval = c(passes = Passes_identifiability, max_delta = max_deltap, max_delta_parm = max_delta_parmp)
            names(rval ) = c("passes","max_delta","max_delta_parm")
            return(rval)
        })
    }
}

#' read $TABLE output from simulation (file name out.dat)
#' and do NCA (Cmax, AUCin,AUClst, from 0 to end.time)
#' if check.identifability, do that with delta_parms as criteria
#' @param run_dir home directory 
#' @param this_sample integer
#' @param NumGroups integer how many groups, e.g., 4 for ABBA
#' @param reference_groups list of arrays, which groups are reference  
#' @param test_groups list of arrays, which groups are  test
#' @param NCA_end_time end time for AUClast and AUCinf 
#' @examples
#' getNCA('c:/runmodels',1,4,c(1,2),c(3,4)),72,0.1)
getNCA = function(run_dir, this_sample, NumGroups, reference_groups, test_groups, NCA_end_time) {
    group_NCA_results <- All_NCA_results <- data.frame(ID = as.integer(), treatment = as.integer(), period = as.integer(), sequence = as.integer(),
        Cmax = as.numeric(), AUCinf = as.numeric(), AUClast = as.numeric())
    data <- read.table(file.path(run_dir, paste0("sim", this_sample), "out.dat"), skip = 1, header = TRUE)
    data <- data %>%
        dplyr::filter(EVID == 0) %>%
        dplyr::filter(DV > 0)
    
    output_file <- file.path(run_dir, paste0("sim", this_sample), paste0("NCAresults", this_sample, ".csv"))
    if (file.exists(output_file)) file.remove(output_file)
  
    #con <- file(output_file, "w") # need to do this to flush file, when using future
    for (this_group in 1:NumGroups) {
        group_data <- data %>%
            dplyr::filter(GROUP == this_group)
        # keep period for each subject, for this group
        period_seq <- group_data %>%
            dplyr::group_by(ID) %>%
            dplyr::distinct(ID, .keep_all = TRUE) %>%
            dplyr::select(ID, OCC, SEQ) %>%
            dplyr::arrange(ID)

        # insert conc=0 at time = 0
        zero_time <- group_data %>%
            dplyr::distinct(ID, .keep_all = TRUE)
        zero_time$TIME <- 0
        zero_time$DV <- 0
        group_data <- rbind(group_data, zero_time) %>%
            dplyr::arrange(ID, TIME)
        conc_obj <- PKNCA::PKNCAconc(group_data, DV ~ TIME | ID)
        data_obj <- PKNCA::PKNCAdata(data.conc = conc_obj, intervals = data.frame(start = 0, end = NCA_end_time, aucinf.obs = TRUE, auclast = TRUE,
            cmax = TRUE))
        results_obj <- PKNCA::pk.nca(data_obj)$result
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
        group_NCA_results <- data.frame(ID = AUCinf$ID, treatment = treatment, period = period_seq$OCC, sequence = period_seq$SEQ, Cmax = CMAX$PPORRES,
            AUCinf = AUCinf$PPORRES, AUClast = AUClast$PPORRES)
         
        All_NCA_results <- rbind(All_NCA_results, group_NCA_results)
    } 
      
    write.csv(All_NCA_results, file = output_file, quote = FALSE, row.names = FALSE)
    # wait for file to be written??, this doesn't seem to help
    count <- 0
    while(!file.exists(output_file) & count < 100){
      Sys.sleep(1)
      count <- count + 1
    }
}



# Right now get_BEQDF(NCAresults1, MetricColumn = 'AUClast', SequenceColumn = 'SEQUENCE') will return AUClast row get_BEQDF(NCAresults1,
# MetricColumn = 'AUCinf', SequenceColumn = 'SEQUENCE') will return AUCinf row sort by subject, treatment and period to get replicate
# number
.get_NCA.Set.ReplicateNumber <- function(NCA.Set) {
    NCA.Set <- data.frame(NCA.Set[order(NCA.Set$ID, NCA.Set$treatment, NCA.Set$period), ])
    NCA.Set$repl <- replicate(nrow(NCA.Set), 1)
    for (i in c(2:nrow(NCA.Set))) {
        if (NCA.Set$ID[i] != NCA.Set$ID[i - 1] | NCA.Set$treatment[i] != NCA.Set$treatment[i - 1]) {
            NCA.Set$repl[i] <- 1
        } else {
            NCA.Set$repl[i] <- NCA.Set$repl[i - 1] + 1
        }
    }

    return(NCA.Set)
}

get_BEQDF <- function(NCA.Set = data.frame(), MetricColumn = "Cmax", SubjectColumn = "ID", TreatmentColumn = "treatment", SequenceColumn = "sequence",
    PeriodColumn = "period", RefValue = "Reference", alpha = 0.05, PartialReplicate = FALSE) {
    RequestedColumnNames <- c(MetricColumn, SubjectColumn, TreatmentColumn, SequenceColumn, PeriodColumn)

    MissingColumns <- is.na(match(RequestedColumnNames, colnames(NCA.Set)))

    if (any(MissingColumns)) {
        stop(paste("Please check the dataset, cannot find the column", RequestedColumnNames[MissingColumns]))
    }

    NCA.Set$logpk <- log(NCA.Set[[MetricColumn]])
    NCA.Set$ID <- NCA.Set[[SubjectColumn]]
    NCA.Set$treatment <- ifelse(NCA.Set[[TreatmentColumn]] == RefValue, "R", "T")
    NCA.Set$period <- NCA.Set[[PeriodColumn]]
    NCA.Set$sequence <- NCA.Set[[SequenceColumn]]

    critbound.s2wR <- get_NCA.Set.RSABE(NCA.Set = NCA.Set, alpha = alpha, PartialReplicate = PartialReplicate)

    ABE <- get_NCA.Set.ABE(NCA.Set, alpha = alpha)
    OverallDF <- cbind.data.frame(MetricColumn = MetricColumn, ABE, critbound.s2wR)
    rownames(OverallDF) <- NULL
    OverallDF$Assessment <- ifelse(OverallDF$swR < 0.294, "Unscaled ABE used", "RSABE used")
    if (OverallDF$swR < 0.294) {
        OverallDF$BE <- ifelse(OverallDF$lower.CL < 0.8 | OverallDF$upper.CL > 1.25, FALSE, TRUE)
    } else {
        OverallDF$BE <- ifelse(OverallDF$pe < 0.8 | OverallDF$pe > 1.25 | OverallDF$critbound > 0, FALSE, TRUE)
    }

    return(OverallDF)
}

# different within variabilities
get_NCA.Set.ABE <- function(NCA.Set, alpha = 0.05) {
    modelABE <- nlme::lme(logpk ~ treatment + period + sequence, subset = !is.na(logpk), random = ~treatment | ID, weights = nlme::varIdent(form = ~treatment),
        data = NCA.Set, method = "REML", na.action = na.exclude, control = list(opt = "optim", msMaxIter = 1000, msMaxEval = 1000))

    # EMA method B modelABE <- lmer(logpk ~ sequence + period + treatment + (1 | ID), data = NCA.Set)

    ConfIntRestults <- confint(emmeans::emmeans(modelABE, pairwise ~ treatment), level = 1 - 2 * alpha)
    Ratio <- exp(-ConfIntRestults$contrasts$estimate) * 100
    # note the difference lower-upper
    lower.CL <- exp(-ConfIntRestults$contrasts$upper.CL)
    upper.CL <- exp(-ConfIntRestults$contrasts$lower.CL)
    return(data.frame(Ratio = Ratio, lower.CL = lower.CL, upper.CL = upper.CL))
}

#' Title
#'
#' @param NCA.Set 
#' @param alpha 
#' @param PartialReplicate 
#'
#' @return
#' @export
#'
#' @examples
get_NCA.Set.RSABE <- function(NCA.Set = data.frame(), alpha = 0.05, PartialReplicate = FALSE) {
    NCA.Set <- .get_NCA.Set.ReplicateNumber(NCA.Set)

    # find not necessary columns
    ColumnNamesToExclude <- colnames(NCA.Set)[!colnames(NCA.Set) %in% c("ID", "sequence", "code2", "logpk")]

    # make a column with treatment and repl, i.e T1,T2, R1,R2
    NCA.Set$code2 <- paste0(NCA.Set$treatment, NCA.Set$repl)
    NCA.Set <- NCA.Set[order(NCA.Set$ID, NCA.Set$period), ]

    # now reshape to wide with cols subj, seq, pk.R1, pk.R2, pk.T1, pk.T2

    NCA.Set_all <- reshape(data = NCA.Set, direction = "wide", v.names = "logpk", idvar = c("ID", "sequence"), timevar = "code2", drop = ColumnNamesToExclude)
    # change logpk.T1, T2 ... to T1, T2 ... to avoid keystrokes
    names(NCA.Set_all) <- sub("logpk.", "", names(NCA.Set_all))

    # now T-R analysis, ilat in the SAS code next will give TR = NA if any T1, T2, R1, R2 is missing must adapted if 3-period TRT|RTR is
    # used
    if (PartialReplicate) {
        NCA.Set_all$TR <- NCA.Set_all$T1 - 0.5 * (NCA.Set_all$R1 + NCA.Set_all$R2)
    } else {
        NCA.Set_all$TR <- 0.5 * (NCA.Set_all$T1 + NCA.Set_all$T2 - NCA.Set_all$R1 - NCA.Set_all$R2)
    }

    # now get rid of subjects not having all 4 periods which have TR = NA
    NCA.Set_ilat <- NCA.Set_all[!is.na(NCA.Set_all$TR), ]  # ilat analysis, ANOVA with seq as effect
    NCA.Set_ilat$sequence <- as.factor(NCA.Set_ilat$sequence)
    # with standard contr.treatment we get not the wanted intercept!  with that the intercept is intercept+seq1
    oc <- options("contrasts")
    on.exit(options(oc))
    options(contrasts = c("contr.sum", "contr.poly"))

    critbound.s2wR <- tryCatch({
        m1 <- lm(TR ~ sequence, data = NCA.Set_ilat)

        # intercept is estimate of ?t-?R
        est <- coef(m1)[1]
        pe <- exp(est)

        # lower, upper 90% CI
        CI <- confint(m1, 1, level = 1 - 2 * alpha)
        dfTR <- m1$df.residual

        # stderr, 'the unknown x' for unbiased estimate of (?T-?R)^2
        stderr <- summary(m1)$coefficients["(Intercept)", "Std. Error"]

        # linearized scABE criterion component x
        x <- est^2 - stderr^2
        boundx <- max(abs(CI))^2

        # now the equivalent of SAS code for R-R, dlat analysis
        NCA.Set_dlat <- NCA.Set_all
        NCA.Set_dlat$RR <- NCA.Set_dlat$R1 - NCA.Set_dlat$R2
        NCA.Set_dlat <- NCA.Set_dlat[!is.na(NCA.Set_dlat$RR), ]
        m2 <- lm(RR ~ sequence, data = NCA.Set_dlat)
        dfRR <- m2$df.residual
        s2wR <- summary(m2)$sigma^2/2

        # progesterone guidance for HVD's
        theta <- ((log(1.25))/0.25)^2

        # linearized scABE criterion component y
        y <- -theta * s2wR
        boundy <- y * dfRR/qchisq(0.95, dfRR)
        swR <- sqrt(s2wR)
        # linearized scABE criterion, 95% upper bound
        crit <- x + y
        critbound <- (x + y) + sqrt(((boundx - x)^2) + ((boundy - y)^2))
        data.frame(swR = swR, pe = pe, critbound = critbound)
    }, error = function(cond) {
        message("Error during linearized scABE criterion calculation:")
        message(cond)
        return(data.frame(swR = NA, pe = NA, critbound = NA))
    })

    return(critbound.s2wR)
}

calc_power = function(run_dir, samp_size){
  
  BEsuccess <- data.frame(Cmax.success = as.logical(), AUCinf.success = as.logical(), AUClast.success = as.logical())
  message(Sys.time(), " Starting statistics for NCA parameters, simulations 1-", samp_size) 
  all_results <- data.frame(Sample = as.integer(),  
                            Cmax_Ratio = as.numeric(), Cmax_lower_CL = as.numeric(), Cmax_upper_CL = as.numeric(), Cmax_swR = as.numeric(), 
                            Cmax_pe = as.numeric(), Cmax_critbound = as.numeric(), Cmax_Assessment = as.numeric(), Cmax_BE = as.logical(),  
                            AUCinf_Ratio = as.numeric(), AUCinf_lower_CL = as.numeric(), AUCinf_upper_CL = as.numeric(), AUCinf_swR = as.numeric(), 
                            AUCinf_pe = as.numeric(), AUCinf_critbound = as.numeric(), AUCinf_Assessment = as.numeric(), AUCinf_BE = as.logical(),  
                            AUClast_Ratio = as.numeric(), AUClast_lower_CL = as.numeric(), AUClast_upper_CL = as.numeric(), AUClast_swR = as.numeric(), 
                            AUClast_pe = as.numeric(), AUClast_critbound = as.numeric(), AUClast_Assessment = as.numeric(), AUClast_BE = as.logical()) 
  
  
  pb <- txtProgressBar(min = 0, max = samp_size, initial = 0)   
  for (this_samp in 1:samp_size) { 
    setTxtProgressBar(pb, this_samp) 
    boolFalse <- FALSE
    count <- 0
    # wait for files to close?? otherwise get Error in file(file, "rt") : cannot open the connection
    while(boolFalse == FALSE & count < 20) {
      data <-  NULL  
      count <- count + 1
      tryCatch({
        data <- read.csv(file.path(run_dir, paste0("sim", this_samp), paste0("NCAresults", this_samp, ".csv")), header = TRUE)
        boolFalse <- TRUE 
      },error = function(e){
        Sys.sleep(1)  
      },finally = {}) # fails will throw error below
    }
    
    Cmax_result <- Cmax_result <- data.frame(MetricColumn = "Cmax", Ratio = -999, lower.CL = -999, upper.CL = -999,  
                                             swR = -999,  pe = -999, critbound = -999, Assessment = -999, BE = -999)
    tryCatch({
      Cmax_result <- get_BEQDF(data %>%  filter(!is.na(Cmax)), MetricColumn = "Cmax", SequenceColumn = "sequence") 
    },error = function(e){
      Cmax_result <- data.frame(MetricColumn = "Cmax", Ratio = -999, lower.CL = -999, upper.CL = -999,  
                                swR = -999,  pe = -999, critbound = -999, Assessment = -999, BE = -999)
    } )
    colnames(Cmax_result) = c("Cmax_MetricColumn","Cmax_Ratio","Cmax_lower_CL","Cmax_upper_.CL","Cmax_swR","Cmax_pe","Cmax_critbound","Cmax_Assessment","Cmax_BE")
    AUClast_result <- AUClast_result <- data.frame(MetricColumn = "AUClast", Ratio = -999, lower.CL = -999, upper.CL=-999,  
                                                   swR = -999,  pe = -999, critbound = -999, Assessment = -999, BE = -999)
    tryCatch({
      AUClast_result <- get_BEQDF(data %>%  filter(!is.na(AUClast)), MetricColumn = "AUClast", SequenceColumn = "sequence") 
    },error = function(e){
      AUClast_result <- data.frame(MetricColumn = "AUClast", Ratio = -999, lower.CL = -999, upper.CL=-999,  
                                   swR = -999,  pe = -999, critbound = -999, Assessment = -999, BE = -999)
    } ) 
    
    colnames(AUClast_result) = c("AUClast_MetricColumn","AUClast_Ratio","AUClast_lower_CL","AUClast_upper_CL","AUClast_swR","AUClast_pe","AUClast_critbound","AUClast_Assessment","AUClast_BE")
    
    AUCinf_result <- AUCinf_result <- data.frame(MetricColumn = "AUCinf", Ratio = -999, lower.CL = -999, upper.CL=-999,  
                                                 swR = -999,  pe = -999, critbound = -999, Assessment = -999, BE = -999)
    tryCatch({
      AUCinf_result <- get_BEQDF(data %>% filter(!is.na(AUCinf)), MetricColumn = "AUCinf", SequenceColumn = "sequence") 
    },error = function(e){
      AUCinf_result <- data.frame(MetricColumn = "AUCinf", Ratio = -999, lower.CL = -999, upper.CL=-999,  
                                  swR = -999,  pe = -999, critbound = -999, Assessment = -999, BE = -999)
    } ) 
    
    colnames(AUCinf_result) <- c("AUCinf_MetricColumn","AUCinf_Ratio","AUCinf_lower_CL","AUCinf_upper_CL","AUCinf_swR","AUCinf_pe","AUCinf_critbound","AUCinf_Assessment","AUCinf_BE")
    all_results <- rbind(all_results,c(Cmax_result,AUCinf_result,AUClast_result))
  
  
  }
close(pb)
return(all_results)
}
make_NCA_plots <- function(BICS, run_dir, samp_size, nmodels, reference_groups, test_groups){
  library(ggplot2)
  BICS = BICS %>% 
    mutate(Samp_num = row_number())
  this_NCAs = NULL
  all_NCAs = data.frame( ID = as.integer(),  treatment = as.character(),	period = as.integer(), 	sequence = as.integer(), 
                         Cmax = as.numeric(),  AUCinf = as.numeric(),	AUClast= as.numeric(),model = as.integer())
  this_model = 1
  for(this_model in 1:nmodels){
    all_NCAs_this_model = data.frame( ID = as.integer(),  treatment = as.character(),	period = as.integer(), 	sequence = as.integer(), 
                                      Cmax = as.numeric(),  AUCinf = as.numeric(),	AUClast= as.numeric(),model = as.integer())
    # gather all Cmax etc from models labeled by model superimpose distributions
    which_models = BICS %>% filter(Best==this_model) %>% select(Samp_num)
    # compile list of NCA with best
    this_samp = which_models$Samp_num[1]
    for(this_samp in which_models$Samp_num){
      this_NCAs =  read.csv(file.path(run_dir,paste0("Sim",this_samp),paste0("NCAresults",this_samp, ".csv")))
      this_NCAs$model = this_model
      all_NCAs_this_model = rbind(all_NCAs_this_model, this_NCAs)
    }
    if(dim(all_NCAs_this_model)[1] > 0){
      all_NCAs = rbind(all_NCAs, all_NCAs_this_model)
      Cmax_plot = ggplot(all_NCAs_this_model,aes(x = Cmax),) + geom_histogram(aes(color=treatment,fill = treatment),position = "identity", alpha=0.4) +
        ggtitle(paste("Cmax, model =",this_model))
      print(Cmax_plot)
      ggsave(file.path(run_dir,paste0("Model_",this_model,"_Cmax_histogram_by_treatment.jpeg")), Cmax_plot, device = "jpeg", width= 9, height= 6)
      AUCinf_plot = ggplot(all_NCAs_this_model,aes(x=AUCinf),) + geom_histogram(aes(color=treatment,fill=treatment),position = "identity",alpha=0.4) +
        ggtitle(paste("AUCinf, model =",this_model))
      print(AUCinf_plot)
      ggsave(file.path(run_dir,paste0("Model_",this_model,"_AUCinf_histogram_by_treatment.jpeg")), AUCinf_plot, device = "jpeg", width= 9, height= 6)
      AUClast_plot = ggplot(all_NCAs_this_model,aes(x=AUClast),) + geom_histogram(aes(color=treatment,fill=treatment),position = "identity",alpha=0.4)+
        ggtitle(paste("AUClast, model =",this_model))
      print(AUClast_plot)
      ggsave(file.path(run_dir,paste0("Model_",this_model,"_AUClast_histogram_by_treatment.jpeg")), AUClast_plot, device = "jpeg", width= 9, height= 6)
    }  
  }
  Cmax_plot = ggplot(all_NCAs,aes(x = Cmax),) + geom_histogram(aes(color=treatment,fill = treatment),position = "identity", alpha=0.4) +
    ggtitle("Cmax, all Models")
  print(Cmax_plot)
  ggsave(file.path(run_dir,"All_models_Cmax_histogram_by_treatment.jpeg"), Cmax_plot, device = "jpeg", width= 9, height= 6)
  AUCinf_plot = ggplot(all_NCAs,aes(x=AUCinf),) + geom_histogram(aes(color=treatment,fill=treatment), position = "identity", alpha=0.4) +
    ggtitle("AUCinf, all Models")
  print(AUCinf_plot)
  ggsave(file.path(run_dir,"All_models_AUCinf_histogram_by_treatment.jpeg"), AUCinf_plot, device = "jpeg", width= 9, height= 6)
  AUClast_plot = ggplot(all_NCAs,aes(x=AUClast),) + geom_histogram(aes(color=treatment,fill=treatment), position = "identity", alpha=0.4) +
    ggtitle("AUClast, all Models")
  print(AUClast_plot)
  ggsave(file.path(run_dir,"All_models_AUClast_histogram_by_treatment.jpeg"), AUClast_plot, device = "jpeg", width= 9, height= 6)
   
  
}

#' run_mbbe_json, reads json file, calls run_mbbe
#' @param Args.json, path to JSON file with arguments
#' @export
#'
#' @examples run_mbbe_json(Args.json)
run_mbbe_json <- function(Args.json) {
  if(!file.exists(Args.json)){
    message(Args.json, " not found, exiting")
  }else{
    Args <- RJSONIO::fromJSON(Args.json) 
    run_mbbe(Args$crash_value, Args$ngroups, Args$reference_groups, Args$test_groups, Args$numParallel, Args$samp_size, Args$run_dir, Args$model_source,
             Args$nmfe_path, Args$delta_parms, Args$use_check_identifiable, Args$NCA_end_time, Args$rndseed, Args$use_simulation_data, Args$simulation_data_path )
  }
}

#' run_mbbe
#'
#' @param crash_value 
#' @param nmodels 
#' @param ngroups 
#' @param reference_groups 
#' @param test_groups 
#' @param numParallel 
#' @param samp_size 
#' @param run_dir 
#' @param model_source 
#' @param nmfe_path 
#' @param delta_parms 
#' @param use_check_identifiable 
#' @param NCA_end_time 
#' @param rndseed 
#' @param use_simulation_data 
#' @param simulation_data_path 
#'
#' @return
#' @export
#'
#' @examples
run_mbbe <- function(crash_value, ngroups, reference_groups, test_groups, numParallel, samp_size, run_dir, model_source, nmfe_path, delta_parms,
    use_check_identifiable, NCA_end_time, rndseed, use_simulation_data, simulation_data_path) {
    message(Sys.time()," Start time\nModel file(s) = ", toString(model_source), "\nreference groups = ", toString(reference_groups), "\ntest groups = ", toString(test_groups))
    message("Bootstrap/Monte Carlo sample size = ", samp_size, "\nnmfe??.bat path = ", nmfe_path, "\nUse_check_identifiability = ", use_check_identifiable) 
    if (use_check_identifiable) {
        message("Delta parameter for use_check_identifiable = ", delta_parms)
    }  
    if (use_simulation_data) {
      message("Simulation data set = ", simulation_data_path)
    } else {
      message("Original analysis data set will be used for simulation")
    }
    
    initial_directory <- getwd()
    path_parents <- split_path(run_dir)
    cur_path <- path_parents[1]
    for (this_parent in 2:length(path_parents)) {
        cur_path <- file.path(cur_path, path_parents[this_parent])
        if (!dir.exists(cur_path)) {
            dir.create(cur_path)
        }
    }
    setwd(run_dir)
    set.seed(rndseed)
     
    msg <- check_requirements(model_source, ngroups, reference_groups, test_groups, nmfe_path, use_check_identifiable, simulation_data_path)
    if (msg$rval) {
        message("Passed requirements check\nCopying source control files from ", toString(model_source), " to ", file.path(run_dir, "modelN"),
            "\n where N is the model number")
        nmodels <- copy_model_files(model_source, run_dir)
        message(Sys.time(), " Sampling data 1-", samp_size," writing data to ", file.path(run_dir, "data_sampM.csv"), " where M is the bootstrap sample number")
        sample_data(run_dir, nmodels, samp_size)
        message(Sys.time(), " Starting bootstrap runs 1-", samp_size," in ", file.path(run_dir, "modelN", "M"), " where N is the model number and M is the sample number",
           "\nProgress bar will appear as models start/complete")
        if (!run.bootstrap(nmfe_path, run_dir, nmodels, samp_size)) {
            message("Failed bootstrap")
        } else {
            # need to wait until all are done, this returns when all are started.
            Sys.sleep(60)  # in case all model are still compiling
            message(Sys.time()," Waiting for bootstrap models to complete")
            wait_for_bs(nmodels, samp_size)
            setwd(run_dir)
            message(Sys.time(), " Getting bootstrap model parameters, samples 1-", samp_size)
            parms <- get_parameters(run_dir, nmodels, samp_size, delta_parms, crash_value, use_check_identifiable)
            base_models <- get_base_model(nmodels)  # get all nmodels base model
            message(Sys.time(), " Constructing simulation  models in ", file.path(run_dir, "SimM"), " where M is the simulation number")
            final_models <- write_sim_controls(run_dir, parms, base_models, samp_size, use_simulation_data, simulation_data_path)  # don't really do anything with final models, already written to disc

            message(Sys.time(), " Running simulation models 1-",samp_size," in ", file.path(run_dir, "SimM"), " where M is the simulation number")
            run_simulations(nmfe_path, run_dir, samp_size)
            Sys.sleep(60)  # in case all model are still compiling, no executable yet
            wait_for_sim(samp_size)
            message(Sys.time()," Calculating NCA parameters for simulations 1-", samp_size, ", writing to ", file.path(run_dir, "SimM", "NCAresultsM"), ",  where M is the simulation number")
            calc_NCA(run_dir, ngroups, reference_groups, test_groups, NCA_end_time, samp_size, numParallel)
            make_NCA_plots(parms$BICS, run_dir, samp_size, nmodels, reference_groups, test_groups)
            # read NCA output and do stats
            all_results <- calc_power(run_dir, samp_size)
            if(!is.null(all_results)){
              write.csv(all_results, file = file.path(run_dir,"All_results.csv"), quote=FALSE)
              Cmax_power <- all_results %>% filter(Cmax_BE != -999) %>%  summarise(Power = mean(Cmax_BE))
              AUClast_power <- all_results %>% filter(AUClast_BE != -999) %>% summarise(Power = mean(AUClast_BE))
              AUCinf_power <- all_results %>% filter(AUCinf_BE != -999) %>% summarise(Power = mean(AUCinf_BE))
              power <- c(Cmax_power, AUCinf_power, AUClast_power)
              print(power)
              write.csv(power,file.path(run_dir,"Power.csv"))
            }else{
              message("Failed in calc_power")
            }
            plan(sequential)
            setwd(initial_directory)
            } 
        }else{ 
          setwd(initial_directory)
          message(msg$msg, " exiting")
    }
}
#
#Args.json <- "u:/fda/mbbe/mbbe/MBBEArgs_mega.json"
Args.json <- "u:/fda/mbbe/mbbe/MBBEArgs_sh.json"
run_mbbe_json(Args.json)
