admm <- function(S, Theta_back, lambda, rho, max_iter = 1000, e_abs = 10^-4, e_rel = 10^-2, quiet = FALSE) {
  p <- nrow(S)
  Theta <- matrix(0, nrow = p, ncol = p)
  Z <- Z_old <- matrix(0, nrow = p, ncol = p)
  U <- matrix(0, nrow = p, ncol = p)

  tolerance <- data.frame()
  is_converged <- FALSE
  score <- 0
  for (i in seq_len(max_iter)) {
    Theta <- update_Theta(S, Z, U, rho)
    Z <- update_Z(Theta, U, Theta_back, lambda, rho)
    U <- update_U(Theta, Z, U)

    n <- p
    e_primal <- n * e_abs + e_rel * max(norm(Theta, "F") , norm(Z, "F"))
    e_dual <- n * e_abs + e_rel * norm(rho * U, "F")
    eval_primal <- norm(Theta - Z, "F")
    eval_dual <- norm(rho * (Z - Z_old), "F")

    df <- data.frame(iter = i, tolerance = c("tolerance", "tolerance", "residual", "residual"),
                     type = c("primal", "dual", "primal", "dual"),
                     value = c(e_primal, e_dual, eval_primal, eval_dual))
    tolerance <- rbind(tolerance, df)
    if (eval_primal <= e_primal && eval_dual <= e_dual) {
      is_converged <- TRUE
      break
    }
    Z_old <- Z
  }
  if (!quiet && !is_converged) warning("rearch max iteration")

  result <- list(Theta = Theta, tolerance = tolerance, is_converged = is_converged)
  class(result) <- "admm"
  result
}

update_Theta <- function(S, Z, U, rho) {
  p <- nrow(S)
  E <- eigen(rho * (Z - U) - S)
  Q <- E$vectors
  lambda <- E$values
  Theta_tilde <- diag((lambda + sqrt(lambda^2 + 4 * rho)) / (2 * rho))
  Theta <- Q %*% Theta_tilde %*% solve(Q)
  Theta
}

update_Z <- function(Theta, U, Theta_back, lambda, rho) {
  soft_thresh_func <- generate_soft_thresh_func(lambda / rho)
  Z <- soft_thresh_func(Theta + U - Theta_back)
  Z
}

generate_soft_thresh_func <- function(lambda) {
  function(x) {
    ifelse(x >= lambda, x - lambda,
           ifelse(x <= -lambda, x + lambda, 0))
  }
}

update_U <- function(Theta, Z, U) {
  U <- Theta - Z + U
  U
}
