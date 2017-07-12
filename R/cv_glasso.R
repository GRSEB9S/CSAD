#' @importFrom glasso glasso
#' @importFrom mvtnorm dmvnorm
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
cv_glasso <- function(df, k=3, lambda=10^seq(-2, 2, length.out = 50), quiet = FALSE) {
  if (!is.data.frame(df)) df <- as.data.frame(df)
  n <- nrow(df)
  Mu <- rep(0, ncol(df))

  partition <- sample(k, n, TRUE)
  inds <- split(seq_len(n), partition)

  if (!quiet) progress_bar <- txtProgressBar(0, length(lambda)+1, style = 3)

  scores <- integer(length(lambda))
  for (i in seq_along(lambda)) {
    rho <- lambda[i]
    logliks <- integer(k)
    for (j in seq_len(k)) {
      df_train <- df[-inds[[j]], ]
      df_test <- df[inds[[j]], ]

      S <- cov(df_train)
      m <- glasso(S, rho=rho)
      loglik <- sum(dmvnorm(df_test, Mu, m$w, log = TRUE))
      logliks[j] <- loglik
    }
    scores[i] <- sum(logliks)
    if (!quiet) setTxtProgressBar(progress_bar, i)
  }

  rho_opt <- lambda[which.max(scores)]

  S <- cov(df)
  model <- glasso(S, rho=rho_opt)
  model$cv <- data.frame(lambda = lambda, loglik = scores)
  model$lambda_opt <- rho_opt
  if (!quiet) setTxtProgressBar(progress_bar, length(lambda)+1)
  model
}
