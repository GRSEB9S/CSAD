#' @export
csad <- function(df_back, df_fore, lambda, rho, max_iter = 1000) {
  # compute background precision matrix
  glasso_back <- cv_glasso(df_back)
  Theta_back <- glasso_back$wi
  # compute empirical covariance matrix for foreground data
  S <- cov(df_fore)
  # search optimum rho
  rho_candidates <- 10^seq(-2, 2, length.out = 10)
  lambda_candidates <- 10^seq(-2, 2, length.out = 10)
  df_loglik <- data.frame()
  for (rho in rho_candidates) {
    message("rho ", rho)
    for (lambda in lambda_candidates) {
      message(lambda)
      loglik <- compute_loglik(df_fore, Theta_back, lambda, rho)
      df_loglik <- rbind(df_loglik, data.frame(rho=rho, lambda=lambda, loglik=loglik))
    }
  }
  ind_opt <- which.max(df_loglik$loglik)
  rho_opt <- df_loglik$rho[ind_opt]
  lambda_opt <- df_loglik$lambda[ind_opt]

  admm <- admm(S, Theta_back, lambda_opt, rho_opt, max_iter)

  result <- list(back = Theta_back, fore = admm$Theta, glasso_back=glasso_back, admm = admm,
                 rho=rho_opt, lambda=lambda_opt, df_loglik=df_loglik)
  class(result) <- "csad"
  result
}

#' @export
bsad <- function(df_back, df_fore) {
  Theta_back <- cv_glasso(df_back)$wi
  Theta_fore <- cv_glasso(df_fore)$wi
  list(back = Theta_back, fore = Theta_fore)
}

compute_loglik <- function(df, Theta_back, lambda, rho) {
  k <- 3
  Mu <- rep(0, ncol(df))

  splitted_df <- split(df, sample(seq_len(k), size = nrow(df), replace = TRUE))
  loglik <- sum(
    sapply(seq_len(k), function(i){
      X_train <- Reduce(rbind, splitted_df[-i])
      X_test <- splitted_df[[i]]
      S <- var(X_train)
      admm <- admm(S, Theta_back, lambda, rho, verbose = FALSE)
      if (!admm$is_converged) return(-Inf)
      Sigma <- solve(admm$Theta)
      loglik <- sum(dmvnorm(X_test, Mu, Sigma, log = TRUE))
      loglik
    })
  )

  loglik
}
