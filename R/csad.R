#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom rBayesianOptimization BayesianOptimization
#' @export
csad <- function(df_back, df_fore,
                 lambda_range = c(0, 0.2), rho_range = c(0, 0.2),
                 measure = c("loglik", "f_measure"), search_method = c("bayesopt", "grid"),
                 bo_init_points = 4, bo_n_iter = 10,
                 max_iter = 1000, quiet = FALSE) {
  # Prepare Arguments -------------------------------------------------------
  if (!is.data.frame(df_back)) df_back <- as.data.frame(df_back)
  if (!is.data.frame(df_fore)) df_fore <- as.data.frame(df_fore)
  measure <- match.arg(measure)
  search_method <- match.arg(search_method)

  # Compute Background Precision Matrix -------------------------------------
  .message("Computing Background Precision Matrix")
  glasso_back <- cv_glasso(df_back, quiet = quiet)
  Theta_back <- glasso_back$wi

  # Search Optimum Params ----------------------------------------------------
  .message("Searching Optimum Hyper Parameters")
  if (measure == "loglik") {
    obj_func <- get_obj_func_loglik(df_fore, Theta_back, k = 3)
  } else {
    obj_func <- get_obj_func_f_measure(df_back, df_fore, k = 3)
  }
  if (search_method == "bayesopt") {
    bo <- BayesianOptimization(
      obj_func,
      bounds = list(lambda = lambda_range, rho = rho_range),
      init_points = bo_init_points, n_iter = bo_n_iter,
      verbose = !quiet
    )
    lambda_opt <- bo$Best_Par["lambda"]
    rho_opt <- bo$Best_Par["rho"]
  } else {
    lambda_opt <- NULL
    rho_opt <- NULL
    max_score <- -Inf
    for (rho in rho_range) {
      for (lambda in lambda_range) {
        score <- obj_func(lambda, rho)$Score
        if (score > max_score) {
          lambda_opt <- lambda
          rho_opt <- rho
          last_score <- score
        }
      }
    }
  }

  # Compute Foreground Precision Matrix -------------------------------------
  S <- cov(df_fore)
  admm <- admm(S, Theta_back, lambda_opt, rho_opt, max_iter)

  result <- list(back = Theta_back, fore = admm$Theta, glasso_back=glasso_back,
                 admm = admm, rho=rho_opt, lambda=lambda_opt)
  class(result) <- "csad"
  result
}

get_obj_func_loglik <- function(df_fore, Theta_back, k = 3) {
  n <- nrow(df_fore)
  Mu <- rep(0, ncol(df_fore))

  partition <- sample(k, n, TRUE)
  inds <- split(seq_len(n), partition)

  min_score <- 0

  function(lambda, rho) {
    logliks <- integer(k)
    for (j in seq_len(k)) {
      df_train <- df_fore[-inds[[j]], ]
      df_test <- df_fore[inds[[j]], ]

      S <- cov(df_train)
      admm <- admm(S, Theta_back, lambda, rho, quiet = TRUE)
      if (!admm$is_converged) return (list(Score = min_score))
      loglik <- compute_loglik(df_fore, admm$Theta)
      logliks[j] <- loglik
    }
    score <- sum(logliks)
    if (score < min_score) min_score <<- score
    list(Score = score)
  }
}

get_obj_func_f_measure <- function(df_back, df_fore, k = 3) {
  n_back <- nrow(df_back)
  n_fore <- nrow(df_fore)

  function(lambda, rho) {
    partition <- sample(k, n_fore, TRUE)
    fore_inds <- split(seq_len(n_fore), partition)
    back_inds <- split(sample(n_back, n_fore), partition)

    predicted <- integer(k * 2)
    for (i in seq_len(k)) {
      df_fore_target <- df_fore[fore_inds[[i]], ]
      df_fore_remain <- df_fore[-fore_inds[[i]], ]
      df_back_target <- df_back[back_inds[[i]], ]
      df_back_remain <- df_back[-back_inds[[i]], ]

      glasso_back <- cv_glasso(df_back_remain, quiet = TRUE)
      Theta_back <- glasso_back$wi

      S <- var(df_fore_remain)
      admm <- admm(S, Theta_back, lambda, rho, quiet = TRUE)
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
}


#' @export
bsad <- function(df_back, df_fore) {
  Theta_back <- cv_glasso(df_back, quiet = TRUE)$wi
  Theta_fore <- cv_glasso(df_fore, quiet = TRUE)$wi
  list(back = Theta_back, fore = Theta_fore)
}


