pcon_desc <- function(dat, decimals = 3) {
  if(dat$design == 'between') {
    if (dat$unbalanced == TRUE) {
      out <- c(dat$ss_b, dat$df_b, dat$ss_c, 1, dat$ss_r, dat$df_r,
             dat$ss_c_unadj, 1)
      dnam <- list(c("between-persons factor", "contrast (weighted)", "residual",
                   "contrast (unweighted)"),
                 c("SS", "df"))
    }
    else {
      out <- c(dat$ss_b, dat$df_b, dat$ss_c, 1, dat$ss_r, dat$df_r)
      dnam <- list(c("between-persons factor", "contrast", "residual"),
                   c("SS", "df"))
    }
    out <- data.frame(matrix(out, ncol = 2, byrow = T, dimnames = dnam))
    out <- format(round(out, decimals), nsmall = decimals)
  }
  else if(dat$design == 'within') {
    out <- c(dat$ss_w, dat$df_w, dat$ss_c, 1, dat$ss_r, dat$df_r)
    dnam <- list(c("within-persons factor", "contrast", "residual"),
                 c("SS", "df"))
    out <- data.frame(matrix(out, ncol = 2, byrow = T, dimnames = dnam))
    out <- format(round(out, decimals), nsmall = decimals)
  }
  else if(dat$design == 'mixed') {
    out <- c(dat$ss_b, dat$df_b, dat$ss_w, dat$df_w, dat$ss_int, dat$df_int,
             dat$ss_c, 1)
    dnam <- list(c("between-persons factor", "within-persons factor",
                   "interaction", "contrast"),
                 c("SS", "df"))
    out <- data.frame(matrix(out, ncol = 2, byrow = T, dimnames = dnam))
    out <- format(round(out, decimals), nsmall = decimals)
  }
  return(out)
}

