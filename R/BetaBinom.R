# dBetaBinom
#' A beta binomial distribution and its helper functions for use in N-mixture
#' extensions
#'
#' \code{dBetaBinom} and \code{dBetaBinom_One} provide a beta binomial
#' distribution nimbleFunction. \code{B} is the beta function. These were written
#' to support the dNmixture beta binomial variations.
#'
#' @name dBetaBinom
#' @aliases B dBetaBinom dBetaBinom_One rBetaBinom rBetaBinom_One
#'
#' @author Ben Goldstein, Lauren Ponisio, and Perry de Valpine
#'
#' @param x vector of integer counts.
#' @param N number of trials, sometimes called "size".
#' @param alpha the first non-negative parameter of the beta distribution for
#'   the binomial probability.
#' @param beta the second non-negative parameter of the beta distribution for
#'   the binomial probability.
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return
#'   probability.
#' @param n number of random draws, each returning a vector of length
#'   \code{len}. Currently only \code{n = 1} is supported, but the argument
#'   exists for standardization of "\code{r}" functions.
#' @param a First non-negative parameter of the beta distribution (called shape1
#'   in base-R)
#' @param b Second non-negative parameter of the beta distribution (called
#'   shape2 in base-R)
#'
#' @details These nimbleFunctions provide distributions that can be
#'     used directly in R or in \code{nimble} hierarchical models (via
#'     \code{\link[nimble]{nimbleCode}} and
#'     \code{\link[nimble]{nimbleModel}}). They were originally written for
#'     the beta binomial N-mixture extensions.
#'
#' The beta binomial distribution is equivalent to a binomial distribution whose
#' probability is itself a beta distributed random variable.
#'
#' The probability mass function of the beta binomial is
#' \code{choose(N, x) * B(x + alpha, N - x + beta) /
#' B(alpha, beta)}.
#'
#' The beta binomial distribution is provided in two forms. \code{dBetaBinom} and
#' \code{rBetaBinom} are used when \code{x} is a vector (i.e. \code{length(x) > 1}),
#' in which case the parameters \code{alpha} and \code{beta} must also be vectors.
#' When \code{x} is scalar, \code{dBetaBinom_One} and \code{rBetaBinom_One} are
#' used.
#'
#' @seealso For beta binomial N-mixture models, see \code{\link{dNmixture}}.
#'   For documentation on the beta function, use \code{?stats::dbeta}
#'
#'
##### Beta binomial support functions #####
#' @rdname dBetaBinom
#' @export
dBetaBinom <- nimbleFunction(
  run = function(x = double(1),
                 N = double(0),
                 alpha = double(1),
                 beta = double(1),
                 log = integer(0, default = 0)) {
    logprob <- 0
    for (i in 1:length(x)) {
      logprob <- logprob +
        B(a = x[i] + alpha[i], b = N - x[i] + beta[i], log = TRUE) -
        B(a = alpha[i], b = beta[ i], log = TRUE) +
        lfactorial(N) - (lfactorial(x[i]) + lfactorial(N - x[i]))
    }

    if (log) return(logprob)
    return(exp(logprob))
    returnType(double(0))
  }
)

#' @rdname dBetaBinom
#' @export
dBetaBinom_One <- nimbleFunction(
  run = function(x = double(0),
                 N = double(0),
                 alpha = double(0),
                 beta = double(0),
                 log = integer(0, default = 0)) {
    logprob <- 0
    logprob <- logprob +
      B(a = x + alpha, b = N - x + beta, log = TRUE) -
      B(a = alpha, b = beta, log = TRUE) +
      lfactorial(N) - (lfactorial(x) + lfactorial(N - x))

    if (log) return(logprob)
    return(exp(logprob))
    returnType(double(0))
  }
)


#' @rdname dBetaBinom
#' @export
rBetaBinom <- nimbleFunction(
  run = function(n = double(0),
                 N = double(0),
                 alpha = double(1),
                 beta = double(1)) {
    p <- numeric(length(alpha))
    for (i in 1:length(alpha)) {
      p[i] <- rbeta(1, alpha[i], beta[i])
    }

    x <- rbinom(length(alpha), N, p)
    return(x)
    returnType(double(1))
  })

#' @rdname dBetaBinom
#' @export
rBetaBinom_One <- nimbleFunction(
  run = function(n = double(0),
                 N = double(0),
                 alpha = double(0),
                 beta = double(0)) {

    p <- rbeta(1, alpha, beta)
    x <- rbinom(1, N, p)
    return(x)
    returnType(double())
  })

#' @rdname dBetaBinom
#' @export
B <- nimbleFunction(
  run = function(a = double(0),
                 b = double(0),
                 log = logical(0)) {
    if (log) return(lgamma(a) + lgamma(b) - lgamma(a + b))
    else return(exp(lgamma(a) + lgamma(b) - lgamma(a + b)))
    returnType(double(0))
  })
