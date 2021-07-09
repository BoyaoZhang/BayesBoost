#' @title cAIC
#' @param y response variable
#' @param y_hat fitted values of the response variable
#' @param X design matrix of fixed effects
#' @param sigma2 variance of epsilon
#' @param G block-diagonal matrix of the var-cov matrix of random effects
#' @param U design matrix of random effects
#' @param q number of random effects variables without intercept
#' @param m number of individuals or clusters
#' @import Matrix
#' @importFrom stats dnorm
#' @export
cAIC <- function(y, y_hat, X, sigma2, G, U, q, m) {
  n <- length(y)

  e <- y - y_hat

  D <- G / sigma2
  D <- forceSymmetric(D)

  Lambdat <- chol(D)
  Lambda <- t(Lambdat)

  Ut <- t(U)

  #-------------------------------------------------------------
  # V0inv
  V0 <- diag(1, nrow = n) + tcrossprod(crossprod(Ut, crossprod(Lambdat)), U)
  # V0 <- diag(1, nrow = n) + U %*% Lambda %*% Lambdat %*% Ut

  # LLt <- tcrossprod(crossprod(Lambda, crossprod(U)), Lambdat) + diag(1, nrow(Lambdat))
  LLt <- Lambdat %*% crossprod(U) %*% Lambda + diag(1, nrow(Lambdat))
  LLt <- forceSymmetric(LLt)
  Lt <- chol(LLt)

  Linv <- solve(t(Lt))

  V0inv <- diag(1, n) - crossprod(Linv %*% (Lambdat %*% Ut))

  #-------------------------------------------------------------
  # A
  XtV0invX <- as.matrix(t(X) %*% V0inv %*% X)
  Rx <- chol(XtV0invX)

  A <- V0inv - crossprod(crossprod(X %*% solve(Rx), V0inv))

  #-------------------------------------------------------------
  # W
  nc <- q + 1
  rowIndices <- rep(1:nc, 1:nc)
  colIndices <- sequence(1:nc)
  LindTemplate <- rowIndices + nc * (colIndices - 1) - choose(colIndices, 2)
  Lind <- rep(LindTemplate, m)

  ind <- Lind
  len <- rep(0, length(Lambda@x))
  Wlist <- list()
  eWelist <- list()
  for (j in 1:length(LindTemplate)) {
    LambdaS <- Lambda
    LambdaSt <- Lambdat
    LambdaS@x <- LambdaSt@x <- len
    LambdaS@x[which(ind == j)] <- LambdaSt@x[which(ind == j)] <- 1
    diagonal <- diag(LambdaS)
    diag(LambdaS) <- diag(LambdaSt) <- 0
    Dj <- LambdaS + LambdaSt
    diag(Dj) <- diagonal
    Wlist[[j]] <- U %*% Dj %*% Ut
    eWelist[[j]] <- as.numeric(e %*% Wlist[[j]] %*% e)
  }

  #----------------------------------------------------------------
  # B
  tye <- as.numeric(t(y) %*% e)

  WAlist <- lapply(Wlist, function(w) {crossprod(t(w), A)})
  tv_WAlist <- lapply(WAlist, function(x) {as.vector(t(x))})
  v_WAlist <- lapply(WAlist, as.vector)

  B <- matrix(0, length(LindTemplate), length(LindTemplate))
  C <- matrix(0, length(LindTemplate), n)

  for (j in 1:length(LindTemplate)) {
    Wj <- Wlist[[j]]
    eWje <- eWelist[[j]]
    C[j, ] <- as.vector(e %*% (Wj %*% A) - eWje * e / (2 * tye))

    for (k in j:length(LindTemplate)) {
      Wk <- Wlist[[k]]
      # WkAWjA <- sum(t(WAlist[[j]]) * WAlist[[k]])
      WkAWjA <- sum(tv_WAlist[[j]] * v_WAlist[[k]])
      eWke <- eWelist[[k]]

      B[j, k] <- B[k, j] <- - tye *
        WkAWjA / (2 * (n - ncol(X))) -
        eWje * eWke / (2 * tye) +
        as.numeric(e %*% Wk %*% (A %*% (Wj %*% e)))
    }
  }

  #-------------------------------------------------------------
  # df
  Lambday <- try(solve(B) %*% C)

  df <- n - sum(diag(A))
  for (j in 1:length(LindTemplate)) {
    df <- df + sum(Lambday[j, ] %*% (A %*% (Wlist[[j]] %*% e)))
  }
  df <- df + 1  # sigma.penalty

  #============================================================
  # cAIC
  bc <- df
  cll <- sum(dnorm(y, y_hat, sd = sqrt(sigma2), log = TRUE))

  cAIC <- -2 * cll + 2 * bc

  return(list(cAIC = cAIC,
              bc = bc,
              cll = cll))
}

