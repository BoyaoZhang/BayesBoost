#' Sample beta from multivariate normal distribution
#' @param y y
#' @param X X
#' @param sigma2 sigma2
#' @param Z design Z matrix
#' @param gamma gamma
#' @keywords internal
sample_beta <- function(y, X, sigma2, Z, gamma) {
  sigma2_beta <- solve((1/sigma2) * (crossprod(X)))
  mu_beta <- tcrossprod(sigma2_beta, t((1/sigma2) * crossprod(X, (y - Z %*% gamma))))
  output <- MASS::mvrnorm(1, mu_beta, sigma2_beta)
  return(output)
}


#' Sample gamma from multivariate normal distribution
#' @param Z Z
#' @param r r
#' @param G G
#' @param y y
#' @param X X
#' @param beta beta
#' @param crossprod_Z crossprod_Z
#' @keywords internal
sample_gamma <- function(Z, r, G, y, X, beta, crossprod_Z) {
  # r <- unique(diag(R))

  solG <- Rfast::spdinv(as.matrix(G))
  # solG <- solve(G)
  Sigma_gamma <- Rfast::spdinv(1/r * crossprod_Z + solG)
  # Sigma_gamma <- solve(tcrossprod(crossprod(Z, solve(R)), t(Z)) + solve(G))
  mu_gamma <- Sigma_gamma %*% (1/r * crossprod(Z, y - X %*% beta))
  # mu_gamma <- tcrossprod(Sigma_gamma, t(1/r * crossprod(Z, y - tcrossprod(X, t(beta)))))
  # mu_gamma <- tcrossprod(Sigma_gamma, t(tcrossprod(crossprod(Z, solve(R)), t(y - tcrossprod(X, t(beta))))))
  mu_gamma <- as.matrix(mu_gamma)

  output <- FastGP::rcpp_rmvnorm(1, Sigma_gamma, mu_gamma)
  output <- as.vector(output)
  # output <- MASS::mvrnorm(1, mu_gamma, Sigma_gamma)
  return(output)
}


#' Sample sigma2 from inverse Gamma distribution
#' @param a a
#' @param b b
#' @param y y
#' @param X X
#' @param beta beta
#' @param Z Z
#' @param gamma gamma
#' @keywords internal
sample_sigma2 <- function(a, b, y, X, beta, Z, gamma) {
  a_tilde <- a + .5 * length(y)  # Regression (Fahrmeir et. al.) p.235

  Xb <- as.vector(X %*% beta)
  Zg <- as.vector(Z %*% gamma)
  cp <- sum((y - Xb - Zg)^2)
  # cp <- crossprod(y - Xb - Zg)
  b_tilde <- b + .5 * cp
  # b_tilde <- as.numeric(b + .5 * cp)
  # b_tilde <- as.numeric(b + .5 * crossprod(y - X %*% beta - Z %*% gamma))
  # b_tilde <- as.numeric(b + .5 * crossprod((y - tcrossprod(X, t(beta)) - tcrossprod(Z, t(gamma))), (y - tcrossprod(X, t(beta)) - tcrossprod(Z, t(gamma)))))
  # output <- MCMCpack::rinvgamma(1, a_tilde, b_tilde)
  output <- 1 / rgamma(1, shape = a_tilde, rate = b_tilde)
  return(output)
}


#' Sample tau2 from inverse Gamma distribution
#' @param a_q a_q
#' @param b_q b_q
#' @param m m
#' @param gamma gamma
#' @keywords internal
sample_tau2 <- function(a_q, b_q, m, gamma) {
  a_q_tilde <- a_q + m / 2
  b_q_tilde <- b_q + .5 * sum(gamma^2)
  # output <- MCMCpack::rinvgamma(1, a_q_tilde, b_q_tilde)
  output <- 1 / rgamma(1, shape = a_q_tilde, rate = b_q_tilde)
  return(output)
}


#' Sample Q from inverse Wishart distribution
#' @param v0 v0
#' @param m m
#' @param gamma gamma
#' @param Lambda0 Lambda0
#' @keywords internal
sample_Q <- function(v0, m, gamma, Lambda0) {
  v <- v0 + m
  l <- NULL
  l <- crossprod(gamma)
  Lambda <- Lambda0 + l
  output <- as.matrix(forceSymmetric(MCMCpack::riwish(v, Lambda)))
  return(output)
}
