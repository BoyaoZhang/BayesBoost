#' @param seed random seed
#' @param n number of replications of each id
#' @param num_id number of individuals/clusters
#' @param beta informative coefficients (including intercept)
#' @param num_base number of baseline variables
#' @param base_mean mean of baseline variables
#' @param base_sd standard deviation of baseline variables
#' @param num_var number of random variables
#' @param var_mean mean of random variables
#' @param var_sd standard deviation of random variables

sim_fun <- function(seed, n, num_id, beta,
                    num_base, base_mean, base_sd,
                    num_var, var_mean, var_sd,
                    tau2, cor_coef, sigma,
                    num_noise,
                    ran_int_model = FALSE) {
  set.seed(seed)
  id.t <- rep(num_id, n)
  id <- rep(1:n, id.t)
  N <- length(id)

  for (i in 1:num_base) {
    nam <- paste0("base", i)
    assign(nam, rep(rnorm(n, base_mean, base_sd), id.t))
  }

  for (i in 1:num_var) {
    nam <- paste0("var", i)
    assign(nam, rnorm(N, var_mean, var_sd))
  }

  X <- as.matrix(as.data.frame(mget(c(paste0("base", 1:num_base), paste0("var", 1:num_var)))))

  ndt <- cor_coef * tau2

  if (ran_int_model) {
    Q <- matrix(tau2)
    colnames(Q) <- rownames(Q) <- c("ranInt")

    Xr <- matrix(1, nrow = nrow(X))
  }
  else {
    Q <- matrix(nrow = 1 + num_var, ncol = 1 + num_var)
    for (i in 1:nrow(Q)) {
      for (j in 1:ncol(Q)) {
        Q[i, j] <- ifelse(i == j, tau2, ndt)
      }
    }
    colnames(Q) <- rownames(Q) <- c("ranInt", paste0("var", 1:num_var))

    Xr <- cbind(1, X[, paste0("var", 1:num_var)])
  }



  Z <- list()
  for (i in 1:n) {
    Z[[i]] <- Xr[id == i, ]
  }
  Z <- as.matrix(Matrix::bdiag(Z))
  gamma <- MASS::mvrnorm(n, rep(0, ncol(Xr)), Q)
  gamma <- as.vector(t(gamma))

  eta <- beta[1] + X %*% beta[-1] + Z %*% gamma
  y <- rnorm(N, eta, sigma)

  if (num_noise != 0) {
    for (i in (1 + num_var):(num_var + num_noise)) {
      assign(paste0("var", i), rnorm(N))
    }

    noise_vars <- as.matrix(as.data.frame(mget(paste0("var", (1+num_var):(num_var+num_noise)))))
    X <- cbind(X, noise_vars)
  }


  # Xr <- X[, pot_ranef]
  # Xnr <- X[, -pot_ranef]
  Xb <- X[, 1:num_base]

  if (!ran_int_model) {
    gamma <- matrix(gamma, ncol = (1 + num_var), byrow = TRUE)
    colnames(gamma) <- c("ranInt", paste0("var", 1:num_var))
  }
  else {
    gamma <- matrix(gamma, ncol = 1)
    colnames(gamma) <- "ranInt"
  }

  return(list(y = y,
              # Xr = Xr,
              # Xnr = Xnr,
              X = X,
              id = id,
              id_constant_vars = Xb,
              gamma = gamma,
              Q = Q))
}

library(BayesBoost)
seed <- 123
data <- sim_fun(seed = seed, n = 50, num_id = 10, beta = c(1, 2, 4, 3, 5),
                num_base = 2, base_mean = 0, base_sd = 1,
                num_var = 2, var_mean = 0, var_sd = 1,
                tau2 = .64, cor_coef = .6, sigma = .4,
                num_noise = 46)
# random intercept
# data <- sim_fun(seed = seed, n = 50, num_id = 10, beta = c(1, 2, 4, 3, 5),
#                 num_base = 2, base_mean = 0, base_sd = 1,
#                 num_var = 2, var_mean = 0, var_sd = 1,
#                 tau2 = .64, cor_coef = 0, sigma = .4,
#                 num_noise = 46, ran_int_model = T)

mod <- bayesboost(data$y, data$X, mstop = 25, id = data$id,
                  id_constant_vars = data$id_constant_vars,
                  t = 30, lambda = .3, pat = 5,
                  max_num_ranef = 5,
                  seed = seed+123)
# mod <- bayesboost(data$y, data$X, mstop = 100, id = data$id,
#                   id_constant_vars = data$id_constant_vars,
#                   t = 30, lambda = .3, pat = 5,
#                   ranef_var = c(3, 4),
#                   seed = seed + 123)
mod$best_mstop




