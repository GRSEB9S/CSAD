#' @importFrom far orthonormalization
#' @export
generate_PSD_matrix <- function(dimension, sparse=0.1) {
  v <- runif(dimension^2, -1, 1)
  v[runif(dimension^2) < 1 - sparse] <- 0
  V <- matrix(v, nrow = dimension, ncol = dimension)
  diag(V) <- runif(dimension, -1, 1)
  O <- orthonormalization(V)
  D <- matrix(0, nrow = dimension, ncol = dimension)
  diag(D) <- rev(sort(runif(dimension, 0, 1)))
  PSD <- t(O) %*% D %*% O
  PSD
}
