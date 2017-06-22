#' @importFrom glasso glasso
#' @importFrom mvtnorm dmvnorm
#' @export
cv_glasso <- function(df, k=3, lambda=10^seq(-10, 10, length.out = 100)) {
  splitted_df <- split(df, sample(seq_len(k), size = nrow(df), replace = TRUE))

  Mu <- rep(0, ncol(df))

  scores <- sapply(lambda, function(rho) {
    sum(
      sapply(seq_len(k), function(i){
        X_train <- Reduce(rbind, splitted_df[-i])
        X_test <- splitted_df[[i]]
        S <- var(X_train)
        m <- glasso(S, rho=rho)
        loglik <- sum(dmvnorm(X_test, Mu, m$w, log = TRUE))
        loglik
      })
    )
  })

  rho_opt <- lambda[which.max(scores)]

  S <- var(df)
  model <- glasso(S, rho=rho_opt)
  model
}
