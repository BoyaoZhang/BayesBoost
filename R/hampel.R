#' @title Hampel Filter
#' @description hampel(x, k = 7, nsigma = 2) applies a Hampel filter to the input vector, x, to detect and remove outliers.
#' For each sample of x, the function computes the median of a window composed of the sample and its six
#' surrounding samples, three per side. It also estimates the standard deviation of each sample about its
#' window median using the median absolute deviation. If a sample differs from the median by more than
#' two standard deviations, it is replaced with the median. If x is a matrix, then hampel treats each
#' column of x as an independent channel.
#' @param x a vector of data points.
#' @param k the window size, default is 7.
#' @param nsigma a number of standard deviations, default is 2.
#' @export

hampel <- function(x, k = 7, nsigma = 2) {
  n <- length(x)
  ind <- NULL
  side_k <- (k - 1) / 2

  S0 <- vector(length = length((1 + side_k):(n - side_k)))
  for (i in (1 + side_k):(n - side_k)) {
    x0 <- median(x[(i - side_k):(i + side_k)])
    S0[i - side_k] <- 1.4826 * median(abs(x[(i - side_k):(i + side_k)] - x0)) # see https://en.wikipedia.org/wiki/Median_absolute_deviation

    if (abs(x[i] - x0) > nsigma * S0[i - side_k]) {
      x[i] <- x0
      ind <- c(ind, i)
    }
  }

  return(list(x = x, ind = ind, S0))
}
