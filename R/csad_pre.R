#' @export
csad_pre <- function(df_back, df_fore, k = 5) {
  rho_candidates <- 10^seq(-2, 2, length.out = 10)
  # rho_candidates <-  0.02782559
  lambda_candidates <- 10^seq(-2, 2, length.out = 10)
  # lambda_candidates <- seq(0, 10, length.out = 20)

  n_back <- nrow(df_back)
  n_fore <- nrow(df_fore)

  result <- data.frame()
  for (rho in rho_candidates) {
    for (lambda in lambda_candidates) {
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

        if (!admm$is_converged) break

        predicted[i * 2 - 1] <-
          which.max(c(compute_loglik(df_back_target, Theta_back),
                      compute_loglik(df_back_target, Theta_fore)))
        predicted[i * 2] <-
          which.max(c(compute_loglik(df_fore_target, Theta_back),
                      compute_loglik(df_fore_target, Theta_fore)))
      }
      if(any(predicted == 0)) break
      actual <- rep(c("back", "fore"), k)
      predicted <- c("back", "fore")[predicted]

      metrics <- as.data.frame(compute_metrics(actual, predicted))
      metrics <- cbind(rho = rho, lambda = lambda, metrics)
      result <- rbind(result, as.data.frame(metrics))

      message("rho: ", rho, "\tlambda: ", lambda,
              "\tprecision: ", metrics$precision,
              "\trecall: ", metrics$recall,
              "\tf_measure: ", metrics$f_measure)
    }
  }
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
