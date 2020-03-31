#' Technical Error of Measurement
#'
#' Compute technical error of measurement and associated normalized versions
#'
#' @export
#' @importFrom stats var na.omit
#'
#' @author David Navega
#'
#' @param x a data.frame of numeric vectors containing multiple replicated
#' observation performed by the same rater/observer or replicated observations
#' performed by different observers.
#' @param digits number of decimal places (rounding)
#'
#' @return a list with the following components:
#' \item{TEM}{Technical error of measurement (in the units of x)}
#' \item{rTEM}{Relative technical error of measurement (as percentage)}
#' \item{R}{Coefficient of Reliability (dimensionless)}
#' \item{n}{an integer with the number of valid observations}
#' \item{k}{an integer indicating the number of observers or replications}
#'
#' @references
#' SJ Ulijaszek, DA Kerr. 1999. Anthropometric measurement error and
#' the assessment of nutritional status. British Journal of Nutrition, 82:
#' 165 - 177.
#'
#' CJ Utermohle, SL Zegura, GM Heathcote. 1983. Multiple Observers, Humidity,
#' and Choice of Precision Statistics: Factors Influencing Craniometric Data
#' Quality. American Journal of Physical Anthropology, 61: 85 -95.
#'
#' @note
#' rTEM, Relative TEM, represents a normalization of TEM using the mean of x.
#' R, Coefficient of Reliability or Variation, represents a normalization of TEM
#' using the variance of x.
#'
#' Missing values are handle by row-wise deletion.
#'
tem <- function(x, digits = 4) {

  if (is.data.frame(x)) {

    m <- na.omit(x)
    n <- nrow(m)
    k <- ncol(m)

    # Technical Error of Measurement, Generalized Formula
    tem <- sqrt(sum(rowSums(m ^ 2) - ((rowSums(m) ^ 2) / k)) / (n * (k - 1)))

    # Relative TEM
    rtem <- tem / mean(colMeans(m)) * 100

    # Coefficient of Reliability or Variation
    r <- (tem ^ 2) / var(unlist(m))

    object <- structure(
      .Data = list(
        TEM = round(tem, digits),
        rTEM = round(rtem, digits),
        R = round(1 - r, digits),
        n = as.integer(n),
        k = as.integer(k)
      ),
      class = "tem"
    )

    return(object)

  } else {
    stop("\n(-) x MUST be a data.frame.")
  }

}


#' Print method for tem
#' @author David Navega
#'
#' @export
#' @noRd
#'
#' @param x an object of class "tem"
#' @param ... ...
#'
print.tem <- function(x, ...) {

  cat("\nSamples:", x$n)
  cat("\nReplicates:", x$k)

  cat("\n\nTechnical Error of Measurement:\n")
  cat("\n TEM:", x$TEM)
  cat("\n Relative TEM (percentage):", x$rTEM)
  cat("\n Coefficient of Reliability:",x$R)

}

