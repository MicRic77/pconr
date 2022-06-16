# calculates inidividual error SS (Rosenthal & Rosnow, 1985)

get_individual_error <- function(dvs, w, subj_n) {
  dvs_temp <- dvs - mean(dvs)
  dvs_temp <- dvs_temp - rowMeans(dvs_temp)
  dvs_temp <- t(t(dvs_temp) - colMeans(dvs_temp))
  l <- rowSums(dvs_temp %*% w)
  ss <- sum(l  * l / sum(w^2))
  df <- subj_n - 1
  ms <- ss / df
  return(c(ss, df, ms))
}

