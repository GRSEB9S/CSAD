#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
csad <- function(df_back, df_fore, max_iter = 1000) {
  # compute background precision matrix
  message("Computing Background Precision Matrix ... ", appendLF = FALSE)
  glasso_back <- cv_glasso(df_back)
  Theta_back <- glasso_back$wi
  message("End")

  # search optimum lambda & rho
  message("Searching Optimum Hyper Parameters ... ")
  rho_candidates <- 10^seq(-2, 2, length.out = 10)
  lambda_candidates <- 10^seq(-2, 2, length.out = 10)
  df_loglik <- data.frame()
  for (rho in rho_candidates) {
    last_loglik <- -Inf
    for (lambda in lambda_candidates) {
      loglik <- compute_loglik(df_fore, Theta_back, lambda, rho)
      df_loglik <- rbind(df_loglik, data.frame(rho=rho, lambda=lambda, loglik=loglik))
      message("rho: ", rho, "\tlambda: ", lambda, "\tloglik: ", loglik)
      if (last_loglik >= loglik) break
      last_loglik <- loglik
    }
  }
  ind_opt <- which.max(df_loglik$loglik)
  rho_opt <- df_loglik$rho[ind_opt]
  lambda_opt <- df_loglik$lambda[ind_opt]
  message("End")

  # compute empirical covariance matrix for foreground data
  S <- cov(df_fore)
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
  if (!is.data.frame(df)) df <- as.data.frame(df)
  Mu <- rep(0, ncol(df))

  splitted_df <- split(df, sample(seq_len(k), size = nrow(df), replace = TRUE))
  loglik <- sum(
    sapply(seq_along(splitted_df), function(i){
      X_train <- Reduce(rbind, splitted_df[-i])
      X_test <- splitted_df[[i]]
      S <- var(X_train)
      admm <- admm(S, Theta_back, lambda, rho, verbose = FALSE)
      if (!admm$is_converged) return(-Inf)
      Sigma <- to_symmetric(solve(admm$Theta))
      loglik <- sum(dmvnorm(X_test, Mu, Sigma, log = TRUE))
      loglik
    })
  )
  loglik
}

to_symmetric <- function(mat) {
  mat[lower.tri(mat)] <- 0
  lower <- t(mat)
  diag(lower) <- 0
  mat + lower
}
