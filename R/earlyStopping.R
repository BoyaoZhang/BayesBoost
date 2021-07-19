#' @title cAIC-based stopping criteria
#' @description Determine the stopping iteration of cAICs based on the existence of
#'   lower values in the following "patience" iterations.
#' @param cAIC a serices of cAICs.
#' @param pat patient paramter, default is 3.
#' @param min_iters the least number of boosting iterations the algorithm should perform,
#'     default is 20. The actual minimum is combined with pat parameter, i.e., "min_iters + pat".
#' @export

earlyStopping <- function(cAIC, pat = 3, min_iters = 20) {
  i <- pat + min_iters
  j <- 0
  v <- Inf
  stop <- i

  while((j < pat) && (i <= length(cAIC) - pat - 1)) {
    i <- i + 1

    if (cAIC[i] < v) {
      j <- 0
      stop <- i
      v <- cAIC[i]
    }
    else {
      j <- j + 1
    }
  }

  return(list(stop = stop))
}

