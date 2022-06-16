pcon_tests <- function(dat, decimals = 3) {
  if(dat$design == 'between') {
    if (dat$unbalanced == TRUE) {
      out <- c(dat$f_c, dat$t_c, NA, dat$p1_c,
             dat$f_r, NA, dat$p2_r, NA,
             dat$f_c_unadj, dat$t_c_unadj, NA, dat$p1_c_unadj)
      dnam <- list(c("contrast (weighted)", "residual",
                     "contrast (unweighted)"),
                 c("F value", "t value", "two-tailed p value",
                   "one-tailed p value"))
    }
    else {
      out <- c(dat$f_c, dat$t_c, NA, dat$p1_c,
               dat$f_r, NA, dat$p2_r, NA)
      dnam <- list(c("contrast", "residual"),
                   c("F value", "t value", "two-tailed p value",
                     "one-tailed p value"))
    }
    out <- data.frame(matrix(out, ncol = 4, byrow = T, dimnames = dnam))
    out <- format(round(out, decimals), nsmall = decimals)
  }
  else if(dat$design == 'within') {
    out <- c(dat$f_c, dat$t_c, NA, dat$p1_c,
             dat$f_c_e_ind, dat$t_c_e_ind, NA, dat$p1_c_e_ind,
             dat$f_r, NA, dat$p2_r, NA)
    dnam <- list(c("contrast (pooled error)", "contrast (individual error)",
                   "residual"),
                 c("F value", "t value", "two-tailed p value",
                   "one-tailed p value"))
    out <- data.frame(matrix(out, ncol = 4, byrow = T, dimnames = dnam))
    out <- format(round(out, decimals), nsmall = decimals)
  }
  else if(dat$design == 'mixed') {
    out <- c(dat$f_c_p, dat$t_c_p, NA, dat$p1_c_p,
             dat$f_c_b, dat$t_c_b, NA, dat$p1_c_b,
             dat$f_c_w, dat$t_c_w, NA, dat$p1_c_w)
    dnam <- list(c("contrast (pooled error)", "contrast (between error)",
                   "contrast (within error)"),
                 c("F value", "t value", "two-tailed p value",
                   "one-tailed p value"))
    out <- data.frame(matrix(out, ncol = 4, byrow = T, dimnames = dnam))
    out <- format(round(out, decimals), nsmall = decimals)
  }
  return(out)
}
