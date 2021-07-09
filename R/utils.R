#' @title Get cluster design matrix Z_i of the random effects
#' @param random_vars random variable matrix without intercept
#' @param n_i number of observations in each cluster
#' @keywords internal
Z_i_fun <- function(random_vars, n_i) {
  # cumsum
  cs <- cumsum(n_i)
  if (!is.null(random_vars)) {
    lapply(1:length(n_i), function(i) {
      if (i == 1)
        Z_i <- random_vars[1:cs[i], ]
      else
        Z_i <- random_vars[(cs[i - 1] + 1):cs[i], ]
      Z_i <- cbind(1, Z_i)
      return(Z_i)
    })
  }
  else {
    lapply(1:length(n_i), function(i) {
      Z_i <- matrix(1, nrow = n_i[i], ncol = 1)
    })
  }
}

#' @title Get cluster design matrix X_i of the fixed effects
#' @param X data matrix without intercept
#' @param n_i number of observations in each cluster
#' @keywords internal
X_i_fun <- function(X, n_i) {
  # cumsum
  cs <- cumsum(n_i)
  lapply(1:length(n_i), function(i) {
    if (i == 1)
      X_i <- X[1:n_i[i], ]
    else
      X_i <- X[(cs[i - 1] + 1):cs[i], ]
    X_i <- cbind(1, X_i)
    return(X_i)
  })
}

#' @title Get clusterwise response y_i
#' @param y response vector
#' @param n_i number of observations in each cluster
#' @keywords internal
y_i_fun <- function(y, n_i) {
  # cumsum
  cs <- cumsum(n_i)
  lapply(1:length(n_i), function(i) {
    if (i == 1)
      y_i <- y[1:n_i[i]]
    else
      y_i <- y[(cs[i - 1] + 1):cs[i]]
    return(y_i)
  })
}


# ## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
# ##   data: a data frame.
# ##   measurevar: the name of a column that contains the variable to be summariezed
# ##   groupvars: a vector containing names of columns that contain grouping variables
# ##   na.rm: a boolean that indicates whether to ignore NA's
# ##   conf.interval: the percent range of the confidence interval (default is 95%)
# summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
#                       conf.interval=.95, .drop=TRUE) {
#
#   # New version of length which can handle NA's: if na.rm==T, don't count them
#   length2 <- function (x, na.rm=FALSE) {
#     if (na.rm) sum(!is.na(x))
#     else       length(x)
#   }
#
#   # This does the summary. For each group's data frame, return a vector with
#   # N, mean, and sd
#   datac <- plyr::ddply(data, groupvars, .drop=.drop,
#                        .fun = function(xx, col) {
#                          c(N    = length2(xx[[col]], na.rm=na.rm),
#                            mean = mean   (xx[[col]], na.rm=na.rm),
#                            sd   = sd     (xx[[col]], na.rm=na.rm)
#                          )
#                        },
#                        measurevar
#   )
#
#   # Rename the "mean" column
#   datac <- rename(datac, c("mean" = measurevar))
#
#   datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
#
#   # Confidence interval multiplier for standard error
#   # Calculate t-statistic for confidence interval:
#   # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
#   ciMult <- qt(conf.interval/2 + .5, datac$N-1)
#   datac$ci <- datac$se * ciMult
#
#   return(datac)
# }

