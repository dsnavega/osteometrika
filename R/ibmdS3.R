#' Information-based Measure of Disagreement
#'
#' Computes Information-based Measure of Disagreement for two or more raters
#'
#'
#' @author David Navega
#' @export
#' @importFrom stats na.omit
#'
#' @param x a data.frame of numeric vectors containing multiple replicated
#' observation performed by the same rater/observer or replicated observations
#' performed by different observers.
#'
#' @param digits number of decimal places (rounding)
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
ibmd <- function(x, digits = 4) {

  if (is.data.frame(x)) {

    m <- na.omit(x)
    n <- nrow(m)
    k <- ncol(m)

    # Compute pairwise IBMD
    ibmd_vector <- vector()
    for(i in 1:(k - 1)) {
      for(j in (i + 1):k) {

        numerator <- abs(x[, c(i, j)][, 1] - x[, c(i, j)][, 2])
        denominator <- apply(x[, c(i, j)], 1, max)
        ibmd_value <- sum(log((numerator / denominator) + 1, base = 2))
        ibmd_vector <- c(ibmd_vector, ibmd_value)

      }
    }

    # Overall IBMD
    ibmd <- sum(ibmd_vector) / (n * ((k * (k - 1)) / 2))

    object <- structure(
      .Data = list(
        IBMD = round(ibmd, digits),
        n = as.integer(n),
        k = as.integer(k)
      ),
      class = "ibmd"
    )

    return(object)

  } else {
    stop("\n(-) x MUST be a data.frame.")
  }

}

#' Print method for ibdm
#' @author David Navega
#'
#' @export
#' @noRd
#'
#' @param x an object of class "ibmd"
#' @param ... ...
#'
print.ibmd <- function(x, ...) {
  cat("\nSamples:", x$n)
  cat("\nReplicates:", x$k)

  cat("\n\nInformation Based Measure of Disagreement:\n")
  cat("\n IBMD:", x$IBMD)
  cat("\n")
}
