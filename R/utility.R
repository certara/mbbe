is.empty <- function(x, mode = NULL, ...)
{
  if (is.null(x)) {
    warning("x is NULL")
    return(FALSE)
  }
  if (is.null(mode))
    mode <- class(x)
  identical(vector(mode, 1), c(x, vector(class(x), 1)))
}

remove_old_files <- function(run_dir, samp_size, model_list){
  msg  <- ""
  # check files that need to be written to, try to delete them
  files_to_remove <-
    file.path(run_dir,
              c(
                "All_results.csv",
                "BICS.csv",
                "Parameters.csv",
                "Power.csv"
              ))
  files_to_remove <-
    c(
      files_to_remove,
      file.path(run_dir, paste0("data_samp", 1:samp_size, ".csv")),
      file.path(
        run_dir,
        paste0("MBBE", 1:samp_size),
        paste0("NCAresults", 1:samp_size, ".csv")
      ),
      file.path(
        run_dir,
        paste0("MBBE", 1:samp_size),
        paste0("MBBE", 1:samp_size, ".mod")
      ),
      file.path(run_dir, paste0("MBBEsim", 1:samp_size),  "OUT.DAT")
    )

  files_to_remove <-
    files_to_remove[file.exists(files_to_remove)]

  count <- 0
  while (any(!file.remove(files_to_remove)) & count < 10) {
    files_to_remove <-
      files_to_remove[file.exists(files_to_remove)]
    count <- count + 1
    sleep(0.1)
  }

  if (length(files_to_remove[file.exists(files_to_remove)]) > 0) {
    msg <- paste("Unable to delete required output file(s) ",
                 paste(files_to_remove, collapse = ", "))
  }

  if (is.null(model_list)) {
    msg <-
      paste(msg,
            "Model list is NULL, error in json file?",
            sep = "\n")
  }else{
    for (this_model in length(model_list)) {
      if (!file.exists(model_list[this_model])) {
        msg <-
          paste(msg,
                paste0("Cannot find ", model_list[this_model]),
                sep = "\n")
        next
      }else{
        control <-
          readLines(model_list[this_model], encoding = "UTF-8", warn = FALSE)

        data_line <- get_block("$DATA", control)
        if (nchar(data_line) == 0) {
          msg <- paste(msg,
                       paste0("In the model file ", model_list[this_model], " DATA block not found."),
                       sep = "\n")
          next
        }

        data_line <-
          gsub("(\\s*\\$DATA\\s*)|(\\s*$)", "", data_line)
        any.quotes <- grep("^\"", data_line)
        if (length(any.quotes) > 0) {
          data_file <-
            sapply(regmatches(
              data_line,
              gregexpr('(\").*?(\")', data_line, perl = TRUE)
            ),
            function(y)
              gsub("^\"|\"$", "", y))[1]

        } else {
          # find first white space
          data_file <- unlist(strsplit(data_line, " "))[1]
        }

        if (!file.exists(data_file)) {
          msg <-
            paste(msg,
                  paste("Cannot find", data_file),
                  sep = "\n")
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
      }
    }
  }
  return(msg)
}
