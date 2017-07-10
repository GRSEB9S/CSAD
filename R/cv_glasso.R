#' @importFrom glasso glasso
#' @importFrom mvtnorm dmvnorm
#' @export
cv_glasso <- function(df, k=3, lambda=10^seq(-2, 2, length.out = 50)) {
  if (!is.data.frame(df)) df <- as.data.frame(df)
  Mu <- rep(0, ncol(df))
  splitted_df <- split(df, sample(k, size = nrow(df), replace = TRUE))
  scores <- sapply(lambda, function(rho) {
    sum(
      sapply(seq_len(k), function(i){
        X_train <- Reduce(rbind, splitted_df[-i])
        X_test <- splitted_df[[i]]
        S <- cov(X_train)
        m <- glasso(S, rho=rho)
        loglik <- sum(dmvnorm(X_test, Mu, m$w, log = TRUE))
        loglik
      })
    )
  })

  rho_opt <- lambda[which.max(scores)]

  S <- cov(df)
  model <- glasso(S, rho=rho_opt)
  model$cv <- data.frame(lambda = lambda, loglik = scores)
  model$lambda_opt <- rho_opt
  model
}
