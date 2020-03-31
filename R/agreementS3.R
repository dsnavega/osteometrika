#' Raw and Chance-corrected Agreement
#'
#' @author David Navega
#' @export
#'
#' @importFrom stats pbinom na.omit
#'
#' @param x an n by data.frame containing n observations performed by m
#' observers for a given attribute of interest.
#' @param weights a character defining the weighting scheme for ordered factors.
#' options are 'linear' or 'quadratic' weights. Default is 'linear'
#' @param random a logical defining if the probability of expected agreement should be
#' computed from totally random allocation model or from the marginals
#' (as in Cohen's Kappa). Default is TRUE (totally random).
#' @param digits number of decimal places (rounding)
#'
#' @return a list with the following components (ordered by relevance):
#' \item{agreement}{a vector with raw and chance-corrected agreement estimates}
#' \item{balanced}{balanced agreement as the mean value of detailed agreement}
#' \item{detailed}{level-specific agreement of a factor}
#' \item{expected}{probability of expected agreement}
#' \item{p.value}{significance probability of a binomial test}
#' \item{n}{number of valid observations/ratings}
#' \item{m}{number of observers/raters}
#' \item{k}{number of levels}
#' \item{probability}{a list of withprobabilities used to compute agreement}
#' \item{ordered}{a logical indicating if ordered factors were analysed}
#' \item{scheme}{weighting scheme applied}
#' \item{weighted}{a logical indicating if analysis is weighted}
#'
#' @details
#' Current implementation generalizes the application of weighting schemes
#' (ordered factors) to raw agreement computation and generalizes
#' chance-corrected agreement to an arbitrary number of raters/rating sessions.
#' Statistical significance is computed from a binomial test against expected
#' probability of agreement and not assuming no agreement (0) as the null
#' hypothesis.
#'

agreement <- function(x, weights = "linear", random = T, digits = 4) {

  # Helpers
  # Cross Tabulate
  cross_tabulate <- function(x) {

    # Exception Handling
    if(!is.data.frame(x))
      stop("\n(-) x must be a data.frame like structure.")

    # Handle tibbles
    if(is.data.frame((x)))
      x <- as.data.frame(x)

    if(!all(sapply(x,is.factor)))
      stop("\n(-) All components of x must be factors.")

    if(any(max(sapply(x, nlevels)) != sapply(x, nlevels)))
      stop("\n(-) Factors contained in x do not have the same number of levels.")

    nraters <- ncol(x)
    counts <- 0
    for(n in 1:(nraters - 1)) {
      for(m in (n + 1):nraters) {
        if(any(levels(x[[n]]) != levels(x[[m]])))
          stop("\n(-) Factor levels for column ", n, " and column ", m," do not match.")
        counts <- counts + table(x[[n]], x[[m]])
      }
    }

    return(counts)

  }

  # Compute Weights for Weighted Joint Probability Computation
  compute_weights <- function(nlevels, type = "linear") {

    levels <- nlevels
    nlevels <- seq_len(levels)
    switch(type,

      linear = {
        weights <- matrix(0, nrow = levels, ncol = levels)
        for (i in nlevels) {
          for(j in nlevels) {
            weights[i, j] <- 1 - (abs(i - j) / (levels - 1))
          }
        }
        weights
      },

      quadratic = {
        weights <- matrix(0, nrow = levels, ncol = levels)
        for (i in nlevels) {
          for(j in nlevels){
            weights[i, j] <- 1 - (abs(i - j) / (levels - 1)) ^ 2
          }
        }
        weights
      },
      {
        stop("\n(-) weighting scheme must be 'linear' or 'quadratic'.")
      }
      )

    return(weights)

  }

  # Algorithm
  scheme <- weights
  # Exception handling
  if(any(sapply(x, class) == "ordered")) {
    if(sum(sapply(x, class) == "ordered") == ncol(x)) {
      ordered <- TRUE
    } else {
      stop("\n[-] x contains a mix of ordered and unordered factors.")
    }
  } else {
    ordered <- FALSE
  }

  # Handle NA by list-wise deletion
  x <- na.omit(x)

  # Compute cross-tabulation to assess agreement
  cross_table <- cross_tabulate(x)
  observed <- cross_table

  # Compute marginals over rows and columns
  p_rows <- rowSums(cross_table)
  p_cols <- colSums(cross_table)

  # Parameters
  nraters <- ncol(x)
  nlevels <- nrow(cross_table)
  nratings <- sum(cross_table)

  # Compute expected agreement for chance-corrected agreement
  expected <- outer(X = p_rows, Y = p_cols, FUN = "*")
  expected <- (expected / nratings ^ 2) * nratings
  if (random) {
    p_random <- 1 / (nlevels * nlevels)
    expected <- matrix(p_random, nrow = nlevels, ncol = nlevels) * nratings
  }

  if (ordered & nlevels > 2) {
    weighted <- TRUE
    # Overall Agreement
    weights <- compute_weights(nlevels = nlevels,type = weights)
    p_observed <- sum(observed  * weights) / nratings
    p_expected <- sum(expected  * weights) / nratings
    p_corrected <- (p_observed - p_expected) / (1 - p_expected)

    # Level-specific Agreement
    omatrix <- (observed * weights) + t(observed * weights)
    omatrix <- omatrix / apply(omatrix, 1, sum)
    ematrix <- (expected * weights) + t(expected * weights)
    ematrix <- ematrix / apply(ematrix, 1, sum)

    ps_observed <- diag(omatrix)
    ps_observed[is.nan(ps_observed)] <- NA

    # Compute Balanced Agreement from level specific agreement
    p_balanced  <- mean(ps_observed, na.rm = T)

  }else{
    weighted <- FALSE
    scheme <- "none"
    # Overall Agreement
    p_observed <- sum(diag(observed)) / nratings
    p_expected <- sum(diag(expected)) / nratings
    p_corrected <- (p_observed - p_expected) / (1 - p_expected)

    # Level-specific Agreement
    omatrix <- observed + t(observed)
    omatrix <- omatrix / apply(omatrix, 1, sum)
    ematrix <- expected + t(expected)
    ematrix <- ematrix / apply(ematrix, 1, sum)
    p_specific <- (omatrix - ematrix) / (1 - ematrix)

    # Level-specific Agreement
    omatrix <- (observed) + t(observed)
    omatrix <- omatrix / apply(omatrix, 1, sum)
    ematrix <- (expected) + t(expected)
    ematrix <- ematrix / apply(ematrix, 1, sum)

    ps_observed <- diag(omatrix)
    ps_observed[is.nan(ps_observed)] <- NA

    # Compute Balanced Agreement from level specific agreement
    p_balanced  <- mean(ps_observed, na.rm = T)

  }

  # Compute statistical significance based on binomial distribution
  p.value <- stats::pbinom(
    q = floor(nratings * p_observed),
    size = nratings,
    prob = p_expected,
    lower.tail = F
  )

  # Prepare function output
  observed <- observed / nratings
  p_marginal <- tabulate(data.matrix(x), nbins = nlevels) / nraters
  p_marginal <- p_marginal / sum(p_marginal)
  names(p_marginal) <- rownames(observed)
  ps_expected <- diag(expected) / nratings
  names(ps_expected) <- rownames(observed)
  agreement <- c(observed = p_observed, corrected = p_corrected)

  structure(
    .Data = list(
      agreement = round(x = agreement, digits = digits),
      balanced = round(x = p_balanced, digits = digits),
      expected = round(x = p_expected, digits = digits),
      p.value = round(p.value, digits = options()$digits),
      detailed = round(ps_observed, digits = digits),
      probability = list(
        joint = round(observed, digits = digits),
        marginal = round(p_marginal, digits = digits),
        expected = round(ps_expected, digits = digits)
      ),
      n = nrow(x),
      m = ncol(x),
      k = nlevels,
      ordered = ordered,
      weighted = weighted,
      scheme = scheme
    ),
    class = "agreement"
  )

}

#' Print method for agreement
#' @author David Navega
#'
#' @export
#' @noRd
#'
#' @param x an object of class "agreement"
#' @param ... ...
#'
print.agreement <- function(x, ...) {

  cat("\nObservations:", x$n)
  cat("\nObservers:", x$m)
  cat("\nLevels:", x$k)
  cat("\nOrdered:", x$ordered)
  cat("\nWeighted:", x$weighted)
  if(x$weighted) {
    cat("\nWeighting:", x$scheme)
  } else {
    cat("\nWeighting:", "none")
  }
  cat("\n\nAgreement:")
  cat("\n   Observed:", x$agreement[1])
  cat("\n  *Balanced:", x$balanced[1])
  cat("\n **Chance-corrected:",x$agreement[2])
  cat("\n   Expected:", x$expected)
  cat("\n   p-value:", x$p.value)
  cat("\n\nLevel-specific Agreement:")
  cat("\n\n   Observed:\n")
  cat("\n")
  print(x$detailed)
  cat("\nLevel Distribution:\n")
  cat("\n")
  print(x$probability$marginal)
  cat("\n(*)  Balanced agreement is the average of level-specific agreements.")
  cat("\n(**) Chance-corrected agreement equates to Cohen's Kappa and its variants.")

}
