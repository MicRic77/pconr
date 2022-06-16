#' Planned contrast analysis
#'
#' The function performs planned contrast analysis using the formulas provides
#' in Maxwell and Delaney (2004) and Rosenthal and Rosnow (1985).
#'
#' @param dat A matrix object. In the case of a between-persons design it needs
#' to include a single between-persons variable in the first column and a
#' dependent variable in the second column. In the case of a
#' @param w A numeric vector object including the contrast weights
#' @param iv_n The number of between-persons variables in data. iv_n can be either 0 or 1.
#'
#' @return A list object including all information association with the planned contrast tests.
#'
#' @references
#' Rosenthal, R., & Rosnow, R. L. (1981). \emph{Contrast analysis: Focused comparisons
#'in the analysis of variance}. Cambridge University Press.
#'
#' Maxwell, S. E., & Delaney, H. D. (2004). \emph{Designing experiments and
#' analyzing data. A model comparison perspective} (2nd ed.). Lawrence Erlbaum.
#'
#' @examples
#' #Example 1 (between-persons design with unequal cell n)
#'
#' cond <- as.factor(c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4))
#' score <- c(3, 5, 1, 7, 3, 2, 2, 9, 4, 5, 7, 4, 9, 4, 1, 8, 5)
#'
#' pcon(cbind(cond, score), w = c(-3, -1, 1, 3), iv_n = 1)
#'
#' @export
#'

pcon <- function(dat, w, iv_n = 0) {

  check_input(dat = dat, iv_n = iv_n)

  out <- list(design = NA, unbalanced = FALSE,
              ss_c = NA, ms_c = NA, f_c = NA, t_c = NA, p2_c = NA, p1_c = NA,
              l_c = NA, psi_c = NA,
              ss_c_unadj = NA, ms_c_unadj = NA, f_c_unadj = NA, t_c_unadj = NA,
              p2_c_unadj = NA, p1_c_unadj = NA, f_c_e_ind = NA, t_c_e_ind = NA,
              p2_c_e_ind = NA, p1_c_e_ind = NA,
              f_c_b = NA, t_c_b = NA, p2_c_b = NA, p1_c_b = NA,
              f_c_b_unvar = NA, t_c_b_unvar = NA, p2_c_b_unvar = NA,
              p1_c_b_unvar = NA,
              f_c_w = NA, t_c_w = NA, p2_c_w = NA, p1_c_w = NA,
              ss_b = NA, df_b = NA, ss_w = NA, df_e = NA,
              ss_int = NA, ss_df = NA,
              ss_e_b = NA, df_e_b = NA, ms_e_b = NA,
              ss_e_b_unvar = NA, df_e_b_unvar = NA, ms_e_b_unvar = NA,
              ss_e_w = NA, df_e_w = NA, ms_e_w = NA,
              ss_e_w_ind = NA, df_e_w_ind = NA, ms_e_w_ind = NA,
              ss_e_m_p = NA, df_e_m_p = NA, ms_e_m_p = NA,
              ss_e_m_ind = NA, df_e_m_ind = NA, ms_e_m_ind = NA,
              ss_r = NA, df_r = NA, ms_r = NA)

  if (iv_n == 0) {
    out$design <- "within"
  }
  else {
    if (dim(dat)[2] == 2) {
      out$design <- "between"
    }
    else {
      out$design <- "mixed"
    }
  }

  if (out$design == "between") {
    iv <- as.factor(dat[, 1])
    subj_n <- table(iv)

    if (min(subj_n) != max(subj_n)) {
      out$unbalanced <- TRUE
    }

    dvs <- dat[, -1]

    w_adj <- table(iv) * (w - (sum(w * table(iv)) / sum(subj_n)))

    c_unadj_res <- get_contrast_var(dvs, w, subj_n, out$design, iv)
    out$ss_c_unadj <- c_unadj_res[2]
    out$psi <- c_unadj_res[1]
    c_res <- get_contrast_var(dvs, w_adj, table(iv), out$design, iv)
    out$ss_c <- c_res[2]
    out$l_c <- c_res[1]

    var_unequal <- get_error_unequalvar(dvs, w, subj_n, iv)
    out$ms_e_b_unvar <- var_unequal[1]
    out$ms_e_b_unvar <- var_unequal[2]

    mse <- get_variances(dvs, iv, subj_n, type = "error_equal_var")[3]
    ssb <- get_between_var(dvs, subj_n, iv)[3]

    out$ms_c_unadj <- out$ss_c_unadj
    out$ms_c <- out$ss_c

    between_aov <- summary(stats::aov(dvs ~ iv))
    out$df_b <- as.numeric(unlist(between_aov)[1])
    out$ss_b <- as.numeric(unlist(between_aov)[3])
    out$df_e_b <- as.numeric(unlist(between_aov)[2])
    out$ss_e_b <- as.numeric(unlist(between_aov)[4])
    out$ms_e_b <- as.numeric(unlist(between_aov)[6])

    out$ss_r <- out$ss_b - out$ss_c
    out$df_r <- out$df_b - 1
    out$ms_r <- out$ss_r / out$df_r

    out$f_c <- out$ms_c / out$ms_e_b
    out$t_c <- sqrt(out$f_c)
    out$p2_c <- 1-stats::pf(out$f_c, 1, out$df_e_b)

    out$f_c_b_unvar <- out$ms_c / out$ms_e_b_unvar
    out$t_c_b_unvar <- sqrt(out$f_c_b_unvar)
    out$p2_c_b_unvar <- 1-stats::pf(out$f_c, 1, out$df_e_b_unvar)

    out$f_c_unadj <- out$ms_c_unadj / out$ms_e_b
    out$t_c_unadj <- sqrt(out$f_c_unadj)

    out$p2_c_unadj <- 1-stats::pf(out$f_c_unadj, 1, out$df_e_b)

    if (out$l_c < 0) {
      out$p1_c_unadj <- 1 - (out$p2_c_unadj / 2)
      out$t_c_unadj <- -1 * out$t_c_unadj
      out$p1_c <- 1 - (out$p2_c / 2)
      out$t_c <- -1 * out$t_c
      out$p1_c_b_unvar <- -1
      out$p1_c_b_unvar <- 1 - (out$p2_c_b_unvar / 2)
    }
    else {
      out$p1_c_unadj <- out$p2_c_unadj / 2
      out$p1_c <- out$p2_c / 2
      out$p1_c_b_unvar <- 1 - (out$p2_c_b_unvar / 2)
    }

    out$f_r <- (out$ss_r / out$df_r) / out$ms_e_b
    out$p2_r <- 1-stats::pf(out$f_r, out$df_r, out$df_e_b)
  }
  else if (out$design == "within") {

    dvs <- dat
    subj_n <- length(dvs[, 1])
    dvs_n <- length(dvs[1,])

    c_res <- get_contrast_var(dvs, w, subj_n, out$design)
    out$ss_c <- c_res[2]
    out$ms_c <- c_res[4]
    out$l_c <- c_res[1]

    dvs_levels <- factor(rep((paste("dv level",1:dvs_n, sep=" ")),
                             rep(subj_n, dvs_n)))
    subj <- factor(rep(paste("subject", 1:subj_n, sep=" ")))
    dvs_vec <- matrix(dvs, ncol=1)
    d <- data.frame(dvs = dvs_vec, subj = subj, time = dvs_levels)

    # ANOVA to determine group.SS and pooled.error.SS
    within_aov <- stats::aov(dvs ~ time + Error(subj / time), data = d)
    out$df_w <- as.numeric(unlist(summary(within_aov))[6])
    out$ss_w <- as.numeric(unlist(summary(within_aov))[8])
    out$df_e_w <- as.numeric(unlist(summary(within_aov))[7])
    out$ss_e_w <- as.numeric(unlist(summary(within_aov))[9])
    out$ms_e_w <- as.numeric(unlist(summary(within_aov))[11])

    ie <- get_individual_error(dvs, w, subj_n)
    out$ss_e_w_ind <- ie[1]
    out$df_e_w_ind <- ie[2]
    out$ms_e_w_ind <- ie[3]

    out$ss_r <- out$ss_w - out$ss_c
    out$df_r <- out$df_w - 1
    out$ms_r <- out$ss_r / out$df_r

    out$f_c <- out$ms_c / out$ms_e_w

    out$t_c <- sqrt(out$f_c)
    out$p2_c <- 1-stats::pf(out$f_c, 1, out$df_e_w)

    out$f_c_e_ind <- out$ms_c / out$ms_e_w_ind
    out$t_c_e_ind <- sqrt(out$f_c_e_ind)
    out$p2_c_e_ind <- 1-stats::pf(out$f_c_e_ind, 1, out$df_e_w_ind)

    if (out$l_c < 0) {
      out$p1_c <- 1 - (out$p2_c / 2)
      out$t_c <- -1 * out$t_c
      out$p1_c_e_ind <- 1 - (out$p2_c_e_ind / 2)
      out$t_c_e_ind <- -1 * out$t_c_e_ind
    }
    else {
      out$p1_c <- out$p2_c / 2
      out$p1_c_e_ind <- out$p2_c_e_ind / 2
    }

    out$f_r <- (out$ss_r / out$df_r) / out$ms_e_w
    out$p2_r <- 1-stats::pf(out$f_r, out$df_r, out$df_e_w)
  }
  else if (out$design == "mixed") {

    iv <- as.factor(dat[, 1])
    iv_nlev <- nlevels(iv)
    dvs <- dat[, -1]
    dvs_nlev <- length(dvs[1,])
    subj_n <- rep(table(iv), each = dvs_nlev)
    subj_ntot <- length(dvs[, 1])

    cell_means <- matrix(rep(NA, iv_nlev * dvs_nlev), nrow = iv_nlev)
    for(i in 1:dvs_nlev) {
      cell_means[, i] <- tapply(dvs[, i], iv, mean)
    }
    cell_means <- as.vector(t(cell_means))

        c_res <- get_contrast_var(cell_means, w, subj_n, out$design)
    out$ss_c <- c_res[2]
    out$ms_c <- c_res[4]
    out$l_c <- c_res[1]

    dvs_level_names <- factor(rep(paste("dv l",1:dvs_nlev, sep=" "), each = subj_ntot))
    subj_names <- factor(rep(paste("subject", 1:subj_ntot, sep=" "), rep = dvs_nlev))
    dvs_vec <- matrix(dvs, ncol=1)
    iv_vec <- rep(iv, dvs_nlev)
    d <- data.frame(iv = iv_vec, dvs = dvs_vec, subj = subj_names, time = dvs_level_names)

    # ANOVA to determine factor.SS and pooled.error.SS
    mixed_aov <- stats::aov(dvs ~ (iv*time) + Error(subj/time), data=d)
    out$df_b <- unlist(summary(mixed_aov))[1]
    out$ss_b <- unlist(summary(mixed_aov))[3]
    out$df_e_b <- unlist(summary(mixed_aov))[2]
    out$ss_e_b <- unlist(summary(mixed_aov))[4]
    out$df_w <- unlist(summary(mixed_aov))[11]
    out$ss_w <- unlist(summary(mixed_aov))[14]
    out$df_int <- unlist(summary(mixed_aov))[12]
    out$ss_int <- unlist(summary(mixed_aov))[15]
    out$df_e_w <- unlist(summary(mixed_aov))[13]
    out$ss_e_w <- unlist(summary(mixed_aov))[16]

    out$ss_e_m_p <- out$ss_e_b + out$ss_e_w
    out$df_e_m_p <- out$df_e_b + out$df_e_w
    out$ms_e_m_p <- out$ss_e_m_p / out$df_e_m_p

    out$ms_e_b <- out$ss_e_b / out$df_e_b
    out$f_c_b <- out$ms_c / out$ms_e_b
    out$t_c_b <- sqrt(out$f_c_b)
    out$p2_c_b <- 1-stats::pf(out$f_c_b, 1, out$df_e_b)

    out$ms_e_w <- out$ss_e_w / out$df_e_w
    out$f_c_w <- out$ms_c / out$ms_e_w
    out$t_c_w <- sqrt(out$f_c_w)
    out$p2_c_w <- 1-stats::pf(out$f_c_w, 1, out$df_e_w)

    out$f_c_p <- out$ms_c / out$ms_e_m_p
    out$t_c_p <- sqrt(out$f_c_p)
    out$p2_c_p <- 1-stats::pf(out$f_c_p, 1, out$df_e_m_p)

    if (out$l_c < 0) {
      out$p1_c_p <- 1 - (out$p2_c_p / 2)
      out$t_c_p <- -1 * out$t_c_p
      out$p1_c_b <- 1 - (out$p2_c_b / 2)
      out$t_c_b <- -1 * out$t_c_b
      out$p1_c_w <- 1 - (out$p2_c_w / 2)
      out$t_c_w <- -1 * out$t_c_w
    }
    else {
      out$p1_c_p <- out$p2_c_p / 2
      out$p1_c_b <- out$p2_c_b / 2
      out$p1_c_w <- out$p2_c_w / 2
    }

  }

  res1 <- pcon_tests(out)
  res2 <- pcon_desc(out)
  res <- list(tests = res1, descriptives = res2)

  return(res)
}
