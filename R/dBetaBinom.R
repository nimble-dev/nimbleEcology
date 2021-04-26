# dBetaBinom
#' A beta binomial distribution and beta function for use in \code{nimble} models
#'
#' \code{dBetaBinom} and \code{dBetaBinom_One} provide a beta binomial
#' distribution that can be used directly from R or in \code{nimble}
#' models. These are also used by beta binomial variations of dNmixture distributions.
#' \code{nimBetaFun} is the beta function.
#'
#' @name dBetaBinom
#' @aliases nimBetaFun dBetaBinom dBetaBinom_One rBetaBinom rBetaBinom_One
#'
#' @author Ben Goldstein and Perry de Valpine
#'
#' @param x vector of integer counts.
#' @param N number of trials, sometimes called "size".
#' @param shape1 shape1 parameter of the beta-binomial distribution.
#' @param shape2 shape2 parameter of the beta-binomial distribution.
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return
#'   probability.
#' @param n number of random draws, each returning a vector of length
#'   \code{len}. Currently only \code{n = 1} is supported, but the argument
#'   exists for standardization of "\code{r}" functions.
#' @param a shape1 argument of the beta function nimBetaFun.
#' @param b shape2 argument of the beta function nimBetaFun.
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
#' \code{choose(N, x) * B(x + shape1, N - x + shape2) /
#' B(shape1, shape2)}, where \code{B(shape1, shape2)} is the beta function.
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
#' @examples
#' # Calculate a beta binomial probability
#' dBetaBinom(x = c(4, 0, 0, 3), N = 10,
#'   shape1 = c(0.5, 0.5, 0.3, 0.5), shape2 = c(0.2, 0.4, 1, 1.2))
#' # Same for case with one observation
#' dBetaBinom_One(x = 3, N = 10, shape1 = 0.5, shape2 = 0.5, log = T)
#' @export

NULL
nimbleOptions(checkNimbleFunction = FALSE)

##### Beta binomial support functions #####
#' @rdname dBetaBinom
#' @export
nimBetaFun <- nimbleFunction(
  run = function(a = double(0),
                 b = double(0),
                 log = logical(0)) {
    if (log) return(lgamma(a) + lgamma(b) - lgamma(a + b))
    else return(exp(lgamma(a) + lgamma(b) - lgamma(a + b)))
    returnType(double(0))
  })

#' @rdname dBetaBinom
#' @export
dBetaBinom <- nimbleFunction(
  run = function(x = double(1),
                 N = double(0),
                 shape1 = double(1),
                 shape2 = double(1),
                 log = integer(0, default = 0)) {
    logprob <- 0
    for (i in 1:length(x)) {
      logprob <- logprob +
        nimBetaFun(a = x[i] + shape1[i], b = N - x[i] + shape2[i], log = TRUE) -
        nimBetaFun(a = shape1[i], b = shape2[ i], log = TRUE) +
        lgamma(N+1) - (lgamma(x[i] + 1) + lgamma(N - x[i] + 1))
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
                 shape1 = double(0),
                 shape2 = double(0),
                 log = integer(0, default = 0)) {
    logprob <- 0
    logprob <- logprob +
      nimBetaFun(a = x + shape1, b = N - x + shape2, log = TRUE) -
      nimBetaFun(a = shape1, b = shape2, log = TRUE) +
      lgamma(N+1) - (lgamma(x+1) + lgamma(N - x + 1))

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
                 shape1 = double(1),
                 shape2 = double(1)) {
    p <- numeric(length(shape1))
    for (i in 1:length(shape1)) {
      p[i] <- rbeta(1, shape1[i], shape2[i])
    }
    x <- rbinom(length(shape1), N, p)
    return(x)
    returnType(double(1))
  })

#' @rdname dBetaBinom
#' @export
rBetaBinom_One <- nimbleFunction(
  run = function(n = double(0),
                 N = double(0),
                 shape1 = double(0),
                 shape2 = double(0)) {
    p <- rbeta(1, shape1, shape2)
    x <- rbinom(1, N, p)
    return(x)
    returnType(double())
  })

nimbleOptions(checkNimbleFunction = TRUE)

