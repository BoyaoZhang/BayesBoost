#' @title BayesBoost algoritm for linear mixed model
#' @param y response variable
#' @param x explanatory variables, it contains all variables (including cluster-constant ones).
#' @param mstop maximal stopping iteration.
#' @param id vector of individuals or clusters.
#' @param id_constant_vars cluster-constant variables.
#' @param t number of MCMC samplings.
#' @param lambda step length or learning rate, default is 0.1.
#' @param hampel_k window-size of Hampel-filter, default is 7, see ?hampel.
#' @param hampel_nsigma number of standard devisions, default is 2, see ?hampel.
#' @param pat "patience" parameter in early stopping, default is 3, see ?earlyStopping.
#' @param max_num_ranef Maximal number of random effect that can be included in the final model,
#' if pot_ranef_var is not NULL, then its default values is the length of pot_ranef_var, else is 2.
#' if max_num_ranef = 0, then execute random intercept model.
#' @param pot_ranef_var A sequence of indexes indicates the potential variables that can be modeled as random effects.
#' @param ranef_var predefined random effects structures, subset of x.
#' @param seed random.seed for reproduction.
#' @importFrom methods new
#' @importFrom stats aggregate density median rgamma
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Matrix bdiag nearPD forceSymmetric
#' @export

bayesboost <- function(y, x, mstop, id, id_constant_vars = NULL,
                       t, lambda = .1, hampel_k = 7, hampel_nsigma = 2, pat = 3,
                       max_num_ranef = ifelse(!is.null(ranef_var), length(ranef_var),
                                              ifelse(is.null(pot_ranef_var), 2, length(pot_ranef_var))),
                       pot_ranef_var = NULL, ranef_var = NULL, seed = NULL) {
  n <- length(y)
  m <- length(unique(id))
  n_i <- table(id)
  n_i_label <- unique(id)

  if (!is.null(ranef_var))
    q <- q_tmp <- length(ranef_var)
  else
    q <- q_tmp <- 0

  if (is.null(colnames(x)))
    colnames(x) <- paste0("x", 1:ncol(x))

  aggr_y <- aggregate(y, by = list(id), mean, na.rm = TRUE)
  y_hat <- matrix(nrow = n, ncol = mstop)
  y_hat[, 1] <- rep(aggr_y$x, n_i)

  x_names <- colnames(x)
  beta <- c(mean(y), rep(0, ncol(x)))
  names(beta) <- c("(Intercept)", x_names)
  beta_mat <- matrix(nrow = ncol(x) + 1, ncol = mstop)
  rownames(beta_mat) <- c("(Intercept)", x_names)
  gamma_sel <- vector(length = mstop)

  mse <- vector(length = mstop)
  best_x <- vector(length = mstop)
  cAIC <- vector(length = mstop)
  bc <- cll <- vector(length = mstop)

  # hat matrix
  hat_mat_list <- lapply(1:ncol(x), function (i) {
    design_mat <- cbind(1, x[, i])
    # (XtX)^(-1)X
    solve(crossprod(design_mat)) %*% t(design_mat)
  })

  g <- gamma <- gamma_list <- list()
  if (!is.null(ranef_var))
    gamma_names <- gamma_names_tmp <- c("ranInt", colnames(x)[ranef_var])
  else
    gamma_names <- gamma_names_tmp <- "ranInt"
  sigma2 <- vector(length = t)
  sigma2_mat <- matrix(nrow = t, ncol = mstop)
  si2 <- r <- 1
  a <- b <- .001
  a_q <- b_q <- .001
  v0 <- 1 + q
  Lambda0 <- Lambda0_tmp <- diag(1, nrow = 1 + q)

  Q <- Q_list <- list()

  Q_mode <- list()
  Q_mode[[1]] <- diag(1, nrow = 1 + q, ncol = 1 + q)
  gamma_mode <- list()
  sigma2_mode <- vector(length = mstop)


  X <- cbind(1, x)

  Z_i <- Z_i_fun(x, n_i)
  Z_all <- as.matrix(bdiag(Z_i))
  #----------------------------
  # orthogonal Z matrix
  # random intercept
  int_column <- seq(1, ncol(Z_all), by = 1 + ncol(x))
  ind_mat <- int_column
  if (!is.null(id_constant_vars)) {
    x_base_int <- cbind(1, id_constant_vars)
    # X(XtX)^(-1)Xt
    xx_int <- x_base_int %*% solve(crossprod(x_base_int)) %*% t(x_base_int)

    Z_int <- Z_all[, int_column]
    Z_int_tilde <- xx_int %*% Z_int
    Z_all[, int_column] <- Z_int - Z_int_tilde
  }


  # random effect
  ranefZ <- lapply(1:ncol(x), function(j) {
    x_ <- matrix(x[, j], ncol = 1)
    xx <- x_ %*% solve(crossprod(x_)) %*% t(x_)
    j_column <- seq(j + 1, ncol(Z_all), by = 1 + ncol(x))
    Z_j <- Z_all[, j_column]
    op <- Z_j - xx %*% Z_j
    return(list(j_column = j_column,
                op = op))
  })
  # for (i in 1:length(ranefZ)) {
  #   Z_all[, ranefZ[[i]]$j_column] <- ranefZ[[i]]$op
  # }
  if (!is.null(ranef_var)) {
    for (i in ranef_var) {
      j_column <- ranefZ[[i]]$j_column
      ind_mat <- rbind(ind_mat, j_column)
      Z_all[, j_column] <- ranefZ[[i]]$op
    }
    Z <- Z_tmp <- Z_all[, as.vector(ind_mat)]
    crossprod_Z <- crossprod_Z_tmp <- as.matrix(crossprod(Z))

    G <- bdiag(replicate(m, bdiag(replicate(length(ranef_var) + 1, bdiag(1)))))
  }
  else {
    Z <- Z_tmp <- Z_all[, int_column]
    crossprod_Z <- crossprod_Z_tmp <- as.matrix(crossprod(Z))

    G <- bdiag(replicate(m, bdiag(1), simplify = FALSE))
  }


  #---------------------------

  pb <- txtProgressBar(0, mstop, style = 3)

  # boosting
  for (s in 1:mstop) {
    # progress bar
    setTxtProgressBar(pb, s)

    # residuals
    if (s == 1)
      res <- y - y_hat[, 1]
    else
      res <- y - y_hat[, s - 1]

    # fit with componentwise least square
    cls <- cls(res, x, hat_mat_list)

    # extract best fitting variable and its coefficients
    best_x[s] <- which.min(cls$mse)
    beta_best_x <- cls$beta_mat[, best_x[s]]

    beta[1] <- beta[1] + lambda * beta_best_x[1]
    beta[1 + best_x[s]] <- beta[1 + best_x[s]] + lambda * beta_best_x[2]

    beta_mat[, s] <- beta

    if (!(best_x[s] %in% match(colnames(id_constant_vars), colnames(x)))) {
      if (is.null(ranef_var)) {
        if (length(gamma_names) - 1 < max_num_ranef) {
          if (is.null(pot_ranef_var) || (!is.null(pot_ranef_var) && (best_x[s] %in% pot_ranef_var))) {
            # select the columns of the Z matrix w.r.t. the best fitting variable
            j_column <- seq(best_x[s] + 1, ncol(Z_all), by = 1 + ncol(x))
            # In order to save the computing TIME,
            # if this variable has never been selected before,
            # replace only the corresponding columns of the best fitting variable in Z matrix
            # with the orthogonal values.
            if (!(best_x[s] %in% best_x[-s]))
              Z_all[, j_column] <- ranefZ[[best_x[s]]]$op

            # If the best fitting variable has not been selected as the random effect,
            # include it as the potential random effect,
            # and update the relevant hyperparameters with temporary values.
            if (!(colnames(x)[best_x[s]] %in% gamma_names)) {
              gamma_names_tmp <- c(gamma_names, colnames(x)[best_x[s]])
              q_tmp <- length(gamma_names_tmp) - 1

              # in the first boosting iteration, use the initialized Q_mode
              if (s == 1)
                Q_mode[[s]] <- bdiag(Q_mode[[s]], 1)
              # otherwise use the Q_mode of the previous iteration
              else
                Q_mode[[s]] <- bdiag(Q_mode[[s - 1]], 1)
              # create potential covariance matrix G with Q_mode
              G <- bdiag(replicate(m, Q_mode[[s]], simplify = FALSE))
              # create potential Lambda0 in inverse Wishart distribution (sample Q)
              Lambda0_tmp <- diag(1, nrow = 1 + q_tmp)

              # the temporary Z matrix is selected from the Z_all matrix
              Z_tmp <- Z_all[, as.vector(rbind(ind_mat, j_column))]
            }
            # else select the temporary Z matrix direct from the Z_all matrix.
            # Note that the selected columns have already been replaced with the orthogonal values.
            else {
              Z_tmp <- Z_all[, as.vector(ind_mat)]
            }

            # ZtZ matrix pre-calculation in order to save computing time
            crossprod_Z_tmp <- as.matrix(crossprod(Z_tmp))
          }
        }
      }
    }




    #-------------------------------------------------
    # Bayesian
    if (!is.null(seed))
      set.seed(seed)
    for (i in 1:t) {
      # sample gamma (vector form)
      g[[i]] <- sample_gamma(Z_tmp, r, G, y, X, beta, crossprod_Z_tmp)
      # gamma (matrix form)
      gamma[[i]] <- t(matrix(g[[i]], nrow = 1 + q_tmp))
      colnames(gamma[[i]]) <- gamma_names_tmp

      # sample sigma2
      r <- sigma2[i] <- sample_sigma2(a, b, y, X, beta, Z_tmp, g[[i]])

      # sample Q
      Q[[i]] <- sample_Q(v0, m, gamma[[i]], Lambda0_tmp)
      rownames(Q[[i]]) <- colnames(Q[[i]]) <- gamma_names_tmp
      # build block-diagonal matrix G from Q
      G <- bdiag_m(replicate(m, Q[[i]], simplify = FALSE))
    }
    # save gamma list
    gamma_list[[s]] <- gamma
    # save sigma2 data
    sigma2_mat[, s] <- sigma2
    # save Q list
    Q_list[[s]] <- Q

    # mode value as the output
    sigma2_mode[s] <- estimate_mode(sigma2)
    # gamma matrix (mode)
    gm <- matrix(apply(array(unlist(g), c(length(g[[1]]), length(g))), 1, estimate_mode), ncol = 1 + q_tmp, byrow = TRUE)
    colnames(gm) <- gamma_names_tmp
    # Q_mode matrix (mode)
    Q_mode_tmp <- as.matrix(nearPD(apply(array(unlist(Q), c(1 + q_tmp, 1 + q_tmp, t)), 1:2, estimate_mode))$mat)

    #------------------------------------------------
    # if (!is.null(ranef_var)) {
    if (length(gamma_names) - 1 < max_num_ranef) {
      # mse Z * gamma
      # mse_Zg <- vector(length = 1 + q_tmp)
      # for each random effect variable, compute the mse
      # for (i in 0:q_tmp) {
      Z_sub <- Z_tmp[, (0:(ncol(Z_tmp) - 1)) %% (1 + q_tmp) == q_tmp]
      mse_Zg <- mean((res - as.vector(X %*% beta) - as.vector(Z_sub %*% gm[, q_tmp + 1]))^2)
      # names(mse_Zg)[q_tmp + 1] <- gamma_names_tmp[q_tmp + 1]
      # }
      # if the smallest mse of the fixed effect smaller than that of the random effect,
      # then discard the best fitting variable as a random effect
      # reset the temporary values (with the values from the previous boosting iteration)
      if (min(cls$mse) <= mse_Zg) {
        q_tmp <- q
        Z_tmp <- Z
        crossprod_Z_tmp <- crossprod_Z
        Lambda0_tmp <- Lambda0
        gamma_names_tmp <- gamma_names
      }
      # else accept the temporary values the newest one
      else {
        # 0 indicates the best fitting variable is random effect instead of fixed effect
        gamma_sel[s] <- gamma_names_tmp[q_tmp + 1]
        # if the best fitting random effect is selected for the first time
        if (!(gamma_names_tmp[q_tmp + 1] %in% gamma_names)) {
          # append it at the end
          gamma_names <- c(gamma_names, gamma_names_tmp[q_tmp + 1])
          # and accept the temporary values
          q <- q_tmp
          Z <- Z_tmp
          crossprod_Z <- crossprod_Z_tmp
          Lambda0 <- Lambda0_tmp
          ind_mat <- rbind(ind_mat, j_column)
        }
        # # else accept directly the temporary values
        else {
          q_tmp <- q
          Z_tmp <- Z
          crossprod_Z_tmp <- crossprod_Z
          Lambda0_tmp <- Lambda0
          gamma_names_tmp <- gamma_names
        }
      }
    }
    # }


    # subset gamma and Q matrix with only informative variables
    gamma_mode[[s]] <- gm[, gamma_names, drop = FALSE]

    rownames(Q_mode_tmp) <- colnames(Q_mode_tmp) <- colnames(gm)
    Q_mode[[s]] <- Q_mode_tmp[gamma_names, gamma_names, drop = FALSE]

    #--------------------------------------------
    # mode values as the starting values of the Bayesian procedure for the next boosting iteration
    r <- sigma2_mode[s]
    G <- bdiag_m(replicate(m, Q_mode[[s]], simplify = FALSE))

    #--------------------------------------------
    # Other statistics
    # MSE
    X_beta <- as.vector(tcrossprod(X, t(beta)))
    Z_gamma <- as.vector(Z %*% as.vector(t(gamma_mode[[s]])))
    y_hat[, s] <- X_beta + Z_gamma
    mse[s] <- mean((y - y_hat[, s])^2)

    # cAIC
    # caic <- cAIC(y, y_hat[, s], X[, which(beta != 0)], sigma2_mode[s], G, Z, q, m)
    caic <- cAICeigen(y, y_hat[, s], X[, which(beta != 0)], sigma2_mode[s], as.matrix(G), as.matrix(Z), q, m)
    cAIC[s] <- caic$cAIC
    bc[s] <- caic$bc
    cll[s] <- caic$cll

  }

  # earlyStopping
  cAIC_hampel <- hampel(cAIC, k = hampel_k, nsigma = hampel_nsigma)$x
  # best_ms <- earlyStopping(cAIC_hampel, pat, best_x)$stop
  best_ms <- earlyStopping(cAIC_hampel, pat)$stop
  x_names_ <- c("ranef", x_names)
  best_x_label <- x_names_[best_x + 1]

  return(list(beta = beta,
              beta_mat = beta_mat,
              mse = mse,
              best_x = best_x,
              best_x_label = best_x_label,
              gamma_sel = gamma_sel,
              gamma_mode = gamma_mode,
              gamma_list = gamma_list,
              sigma2_mode = sigma2_mode,
              sigma2_mat = sigma2_mat,
              Q_list = Q_list,
              Q_mode = Q_mode,
              y_hat = y_hat,
              cAIC = cAIC,
              cAIC_hampel = cAIC_hampel,
              bc = bc,
              cll = cll,
              best_mstop = best_ms))
}



#' @title Fast version of Matrix:: .bdiag()
#' @param lmat matrix
#' @keywords internal
bdiag_m <- function(lmat) {
  ## Copyright (C) 2016 Martin Maechler, ETH Zurich
  if(!length(lmat)) return(new("dgCMatrix"))
  stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
            (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
            all(vapply(lmat, dim, integer(2)) == k)) # all of them
  N <- length(lmat)
  if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- as.integer(N * k)
  ## result: an   M x M  matrix
  new("dgCMatrix", Dim = c(M,M),
      ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
      i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
      p = k * 0L:M,
      x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
}

#' @title Componentwise least square
#' @param y y
#' @param x x
#' @param hat_mat_list list of hat matrix
#' @keywords internal
cls <- function(y, x, hat_mat_list) {
  mse <- vector(length = ncol(x))
  op <- lapply(1:ncol(x), function (i) {
    beta_hat <- tcrossprod(hat_mat_list[[i]], t(y))
  })
  beta_mat <- matrix(unlist(op), nrow = 2, ncol = length(op))
  pred <- sapply(1:ncol(beta_mat), function (i) {beta_mat[1, i] + beta_mat[2, i] * x[, i]})
  mse <- sapply(1:ncol(pred), function (i) {mean((y - pred[, i])^2)})
  return(list(beta_mat = beta_mat,
              mse = mse))
}

#' @title get posterior mode
#' @param x x
#' @keywords internal
estimate_mode <- function(x) {
  d <- density(x, from = min(x), to = max(x))
  d$x[which.max(d$y)]
}


