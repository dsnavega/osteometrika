#' Kendall's W Coefficient of Concordance
#'
#' @author David Navega
#' @export
#'
#' @importFrom stats cor pchisq na.omit p.adjust
#'
#' @param x an n by m matrix or data.frame containing n observations performed
#' by m observers for a given attribute of interest
#' @param pairwise a logical indicating if pairwise detailed analysis should be
#' performed if more than two observers are present (columns of x.).
#' Default is TRUE
#' @param digits number of decimal places (rounding)
#'
#' @return a list with the following components:
#' \item{n}{number of valid observations/ratings}
#' \item{m}{number of observers/raters}
#' \item{estimate}{the estimate of Kendall's W}
#' \item{statistic}{test statistic for the Chi-Squared distribuition}
#' \item{p.value}{p-value for the computed test}
#' \item{pairwise}{a data.frame with pairwise analysis between observers/raters}
#'
#' @note Current implementations computes Kendall's W from Spearman's Rho
#'
concordance <- function(x, pairwise = T, digits = 4) {

  x <- na.omit(data.matrix(x))
  n <- nrow(x)
  m <- ncol(x)

  # Compute Kendall's W from Spearman's Rho
  rho <- cor(x, method = "spearman")
  pairs <- upper.tri(x = rho, diag = F)
  estimate <- (mean(rho[pairs]) + 1) / m

  # Statistical significance based on Chi-Squared Distribution
  statistic <- m * (n - 1) * estimate
  pvalue <- pchisq(q = statistic, df = n - 1, lower.tail = F)

  # Pairwise Analysis
  pw_matrix <- NULL
  if(pairwise & m > 2) {

    for(i in 1:(m - 1)) {
      for(j in (i + 1):m) {
        estimate <- (rho[i, j] + 1) / 2
        statistic <- 2 * (n - 1) * estimate
        pvalue <- p.adjust(pchisq(q = statistic, df = n - 1, lower.tail = F))
        pw_matrix <- rbind(pw_matrix, c(i, j, estimate, statistic, pvalue))
      }
    }
    pw_matrix <- data.frame(pw_matrix)
    colnames(pw_matrix) <- c("i", "j", "estimate", "statistic", "p.value")

  }

  object <- structure(
    .Data = list(
      n = n,
      m = m,
      estimate = round(estimate, digits),
      statistic = round(statistic, digits),
      p.value = round(pvalue, options()$digits),
      pairwise = round(pw_matrix, digits)
    ),
    class = "concordance"
  )

  return(object)

}



#' Print method for concordance
#' @author David Navega
#'
#' @export
#' @noRd
#'
#' @param x an object of class "concordance"
#' @param ... ...
#'
print.concordance <- function(x, ...) {

  cat("\nSamples:", x$n)
  cat("\nReplicates:", x$m)

  cat("\n\nKendall's W:")
  cat("\n Estimate (W):", x$estimate)
  cat("\n Statistic (Chi-squared):", x$statistic)
  cat("\n p-value:", x$p.value)

  cat("\n\nPairwise Analysis:\n\n")
  pairwise <- x$pairwise
  rownames(pairwise) <- apply(
    pairwise[,1:2], 1, function(x) paste(x,collapse = " - ")
  )

  print(pairwise)

}
