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

get_pid <- function(exeName) {
  tasklist_raw <- system("tasklist", intern = TRUE)
  tasklist_raw_mat <-
    suppressWarnings(do.call(rbind, strsplit(tasklist_raw[-1:-3],
                                              "   +")))
  image_name <- tasklist_raw_mat[, 1]
  tasklist_raw_mat_fix_pids <- do.call(rbind, strsplit(tasklist_raw_mat[,
                                                                        2], " "))
  tasklist <- data.frame(image_name, PID = tasklist_raw_mat_fix_pids, stringsAsFactors = FALSE)
  process_pid <- tasklist[exeName, "PID"]
}
