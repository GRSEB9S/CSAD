#' @export
csad <- function(df_back, df_fore, lambda, rho, max_iter = 1000) {
  Theta_back <- cv_glasso(df_back)$wi
  S <- cov(df_fore)
  admm <- admm(S, Theta_back, lambda, rho, max_iter)
  result <- list(back = Theta_back, fore = admm$Theta, admm = admm)
  class(result) <- "csad"
  result
}

#' @export
bsad <- function(df_back, df_fore) {
  Theta_back <- cv_glasso(df_back)$wi
  Theta_fore <- cv_glasso(df_fore)$wi
  list(back = Theta_back, fore = Theta_fore)
}
