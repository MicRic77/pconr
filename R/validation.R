# function to check the pcon input ----

check_input <- function(dat, iv_n) {

  error <- F
  error_msg <- NULL

  if(!is.matrix(dat)){
    error <- T
    error_msg <- append(error_msg, "Invalid input: dat needs to be a matrix
                        object.")
  }
  else {
    if(dim(dat)[2] < 2) {
      error = T
      error_msg <- append(error_msg, "Invalid input: dat needs to have at least
      two columns.f")
    }
  }

  if (is.vector(iv_n) & length(iv_n) == 1) {
    if (! iv_n %in% c(0, 1)) {
      error = T
      error_msg <- append(error_msg, "Invalid input: iv_n needs to be either 0
(for within-subject designs) or 1 (for between-subject or mixed designs.")
    }
  }
  else {
    error = T
    error_msg <- append(error_msg, "Invalid input: iv_needs to be a single
number")
  }

  if (error == T) {
    stop(error_msg)
  }
}

