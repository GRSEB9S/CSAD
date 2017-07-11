.message <- function(..., domain = NULL, appendLF = TRUE) {
  quiet <- get("quiet", envir = parent.frame())
  if (!quiet) {
    message(..., domain = domain, appendLF = appendLF)
  }
}

#' @export
compute_loglik <- function(df, Theta) {
  Mu <- rep(0, ncol(df))
  Sigma <- to_symmetric(solve(Theta))
  sum(dmvnorm(df, Mu, Sigma, log = TRUE))
}

to_symmetric <- function(mat) {
  mat[lower.tri(mat)] <- 0
  lower <- t(mat)
  diag(lower) <- 0
  mat + lower
}

put_list <- function (list, envir = parent.frame()) {
  message(paste(sprintf("%s: %s", names(list), list),
                collapse = ", "))
}
