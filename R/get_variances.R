get_variances <- function(dvs, iv = NULL, subj_n, w = NULL, type) {
  if (type == "error_equal_var") {
    v <- tapply(dvs, iv, var)
    df <- sum(subj_n) - nlevels(iv)
    ss <- sum((subj_n - 1) * v)
    ms <- ss / df
  }
  else if (type == "error_unequal_var") {
    # according to Maxwell and Delaney
    v <- tapply(dvs, iv, var)
    k <- sum(w^2 / subj_n)
    ms <- sum(k * v) / sum(k)
    k <- w^2 * v / subj_n
    df <- sum(k)^2 / sum(k^2 / (subj_n - 1))
    ss <- ms * df
  }
  else if (type == "error_individual") {
    # according to Rosenthal and Rosnow (1985)
    dvs_temp <- dvs - mean(dvs)
    dvs_temp <- dvs_temp - rowMeans(dvs_temp)
    dvs_temp <- t(t(dvs_temp) - colMeans(dvs_temp))
    psi <- colSums(dvs_temp %*% w)
    ss <- sum(psi  * psi / sum(w^2))
    df <- subj_n - psi
    ms <- ss / df
  }
  else if (type == "between") {
    m <- tapply(dvs, iv, mean)
    df <- nlevels(iv) - 1
    ss <-   sum(subj_n * (m - mean(m))^2)
    ms <- ss_b / df_b
  }

  return(c(ss, df, ms))
}

get_error_var <- function(dvs, subj_n, iv) {

  v <- tapply(dvs, iv, var)
  df_e <- sum(subj_n) - nlevels(iv)
  ss_e <- sum((subj_n - 1) * v)
  ms_e <- ss_e / df_e

  return(c(ss_e, df_e, ms_e))
}

get_error_unequalvar <- function(dvs, w, subj_n, iv) {

  v <- tapply(dvs, iv, var)
  k <- sum(w^2 / subj_n)
  ms_e <- sum(k * v) / sum(k)

  k <- w^2 * v / subj_n
  df_e <- sum(k)^2 / sum(k^2 / (subj_n - 1))

  return(c(ms_e, df_e))
}



get_between_var <- function(dvs, subj_n, iv) {

  m <- tapply(dvs, iv, mean)
  df_b <- nlevels(iv) - 1
  ss_b <-   sum(subj_n * (m - mean(m))^2)
  ms_b <- ss_b / df_b

  return(c(ss_b, df_b, ms_b))
}

get_contrast_var <- function(dvs, w, subj_n, design, iv = NULL) {

  if (design == "between") {
    psi <- sum(tapply(dvs, iv, mean) * w)
  }
  else if (design == "within") {
    psi <- sum(apply(dvs, 2, mean) * w)

  }
  else if (design == "mixed") {
    psi <- sum(dvs * w)
  }

  df_c <- 1
  ss_c <- psi^2  / (sum(w^2 / subj_n))
  ms_c <- ss_c

  return(c(psi, ss_c, df_c, ms_c))
}


