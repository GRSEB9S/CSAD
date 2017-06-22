#' @export
csad <- function(df_back, df_fore) {
  Theta_back <- cv_glasso(df_back)$wi
  S <- cov(df_fore)
  Theta_fore <- admm(S, Theta_back, lambda = .3, rho = .1)
  list(back = Theta_back, fore = Theta_fore)
}

#' @export
bsad <- function(df_back, df_fore) {
  Theta_back <- cv_glasso(df_back)$wi
  Theta_fore <- cv_glasso(df_fore)$wi
  list(back = Theta_back, fore = Theta_fore)
}
