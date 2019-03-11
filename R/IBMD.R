#' Information-based Measure of Disagreement
#'
#' Computes Information-based Measure of Disagreement for two or more raters
#'
#' @export
#' @importFrom stats complete.cases
#'
#' @author David Navega
#'
#' @param x a data.frame of numeric vectors containing multiple replicated
#' observation performed by the same rater/observer or replicated observations
#' performed by different observers.
#'
#' @return a list with the following components:
#' \item{IBMD}{a numeric vector with the value of the IBMD}
#' \item{n}{an integer with the number of valid observations}
#' \item{k}{an integer indicating the number of observers or replications}
#'
#' @references
#' T Henriques, L Antunes, J Bernardes, M Matias, D Sato, C Costa-Santos. 2013.
#' Information-based measure of disagreement for more than two observers:
#' a useful tool to compare the degree of observer disagreement. BMC Medical
#' Research Methodology 13:47. https://doi.org/10.1186/1471-2288-13-47
#'
#' C Costa-Santos, L Antunes, A Souto, J Bernardes. 2010. Assessment of
#' Disagreement: A New Information-Based Approach. Annals of Epidemiology,
#' 20(7): 555 - 561. https://doi.org/10.1016/j.annepidem.2010.02.011
#'
#' @note
#' Missing values are handle by row-wise deletion.
#'
IBMD <- function(x) {

  condition <- is.data.frame(x)
  if(condition) {

    # Intialize ----
    no_na <- complete.cases(x)
    x <- x[no_na, ]
    n <- nrow(x)
    k <- ncol(x)

    # Compute pairwise IBMD ----
    ibmd_vector <- vector()
    for(i in 1:(k - 1)) {
      for(j in (i + 1):k) {

        numerator <- abs(x[, c(i, j)][, 1] - x[, c(i, j)][, 2])
        denominator <- apply(x[, c(i, j)], 1, max)
        ibmd_value <- sum(log((numerator / denominator) + 1, base = 2))
        ibmd_vector <- c(ibmd_vector, ibmd_value)

      }
    }

    # Overall IBMD ----
    ibmd <- sum(ibmd_vector) / (n * ((k * (k - 1)) / 2))

    # return
    rout <- list(
      IBMD = ibmd,
      n = as.integer(n),
      k = as.integer(k)
    )

    return(rout)

  } else {
    stop("[-] x MUST be a data.frame.")
  }

}
