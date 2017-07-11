#' @importFrom rBayesianOptimization BayesianOptimization
#' @export
csad_pre <- function(df_back, df_fore, k = 5, lambda = c(0, 0.2), rho = c(0, 0.2),
                     init_points = 4, n_iter = 10) {
  if(is.data.frame(df_back)) df_back < as.data.frame(df_back)
  if(is.data.frame(df_fore)) df_fore < as.data.frame(df_fore)

  n_back <- nrow(df_back)
  n_fore <- nrow(df_fore)

  obj_func <- function(lambda, rho) {
    partition <- sample(k, n_fore, TRUE)
    fore_inds <- split(seq_len(n_fore), partition)
    back_inds <- split(sample(n_back, n_fore), partition)

    predicted <- integer(k * 2)
    for (i in seq_len(k)) {
      df_fore_target <- df_fore[fore_inds[[i]], ]
      df_fore_remain <- df_fore[-fore_inds[[i]], ]
      df_back_target <- df_back[back_inds[[i]], ]
      df_back_remain <- df_back[-back_inds[[i]], ]

      glasso_back <- cv_glasso(df_back_remain)
      Theta_back <- glasso_back$wi

      S <- var(df_fore_remain)
      admm <- admm(S, Theta_back, lambda, rho, quiet = FALSE)
      Theta_fore <- admm$Theta

      if (!admm$is_converged) return(list(Score=0))

      predicted[i * 2 - 1] <-
        which.max(c(compute_loglik(df_back_target, Theta_back),
                    compute_loglik(df_back_target, Theta_fore)))
      predicted[i * 2] <-
        which.max(c(compute_loglik(df_fore_target, Theta_back),
                    compute_loglik(df_fore_target, Theta_fore)))
    }

    actual <- rep(c("back", "fore"), k)
    predicted <- c("back", "fore")[predicted]

    metrics <- compute_metrics(actual, predicted)
    list(Score = metrics$f_measure + rnorm(1, 0, 1e-3), metrics = metrics)
  }
  bo <- BayesianOptimization(
    obj_func,
    bounds = list(lambda = lambda, rho = rho),
    init_points = init_points, n_iter = n_iter
  )
  rho_best <- bo$Best_Par["rho"]
  lambda_candidates <- seq(lambda[1], lambda[2], length.out = 10)

  df <- data.frame()
  for (lambda in lambda_candidates) {
    metrics <- obj_func(lambda, rho_best)$metrics
    df_metrics <- cbind(rho = rho_best, lambda = lambda, as.data.frame(metrics))
    df <- rbind(df, df_metrics)
    put(as.list(df_metrics))
  }

  result <- list(rho_best = rho_best, metrics = df)
  class(result) <- "csad_pre"
  result
}



compute_metrics <- function(actual, predicted) {
  confusion_matrix <- create_confusion_matrix(actual, predicted)
  true_positive <- confusion_matrix[1, 1]
  false_positive <- confusion_matrix[2, 1]
  false_negative <- confusion_matrix[1, 2]
  precision <- true_positive / (true_positive + false_positive)
  recall <- true_positive / (false_negative + true_positive)
  f_measure <- 2 * precision * recall / (precision + recall)
  list(precision = precision, recall = recall, f_measure = f_measure)
}


create_confusion_matrix <- function(actual, predicted) {
  confusion_matrix <- table(actual, predicted)
  if(ncol(confusion_matrix) == 1) {
    if(colnames(confusion_matrix) == "fore") {
      confusion_matrix <- cbind(back=c(0,0), confusion_matrix)
    } else {
      confusion_matrix <- cbind(confusion_matrix, fore=c(0,0))
    }
  }
  confusion_matrix
}
