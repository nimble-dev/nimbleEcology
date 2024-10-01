# dBetaBinom
#' A beta binomial distribution and beta function for use in \code{nimble} models
#'
#' \code{dBetaBinom_v} and \code{dBetaBinom_s} provide a beta binomial
#' distribution that can be used directly from R or in \code{nimble}
#' models. These are also used by beta binomial variations of dNmixture distributions.
#' \code{nimBetaFun} is the beta function.
#'
#' @name dBetaBinom
#' @aliases nimBetaFun dBetaBinom_v dBetaBinom_s rBetaBinom_v rBetaBinom_s
#'
#' @author Ben Goldstein and Perry de Valpine
#'
#' @param x vector of integer counts.
#' @param N number of trials, sometimes called "size".
#' @param shape1 shape1 parameter of the beta distribution.
#' @param shape2 shape2 parameter of the beta distribution.
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return
#'   probability.
#' @param n number of random draws, each returning a vector of length
#'   \code{len}. Currently only \code{n = 1} is supported, but the argument
#'   exists for standardization of "\code{r}" functions.
#' @param len length of \code{x}.
#' @param a shape1 argument of the beta function.
#' @param b shape2 argument of the beta function.
#'
#' @details These nimbleFunctions provide distributions that can be used
#'   directly in R or in \code{nimble} hierarchical models (via
#'   \code{\link[nimble]{nimbleCode}} and \code{\link[nimble]{nimbleModel}}).
#'   They are used by the beta-binomial variants of the N-mixture distributions
#'   (\code{\link{dNmixture}}).
#'
#' The beta binomial is the marginal distribution of a  binomial distribution whose
#' probability follows a  beta distribution.
#'
#' The probability mass function of the beta binomial is
#' \code{choose(N, x) * B(x + shape1, N - x + shape2) /
#' B(shape1, shape2)}, where \code{B(shape1, shape2)} is the beta function.
#'
#' \code{nimBetaFun(shape1, shape2)} calculates \code{B(shape1, shape2)}.
#'
#' The beta binomial distribution is provided in two forms. \code{dBetaBinom_v} and
#' is when \code{shape1} and \code{shape2} are vectors.
#' \code{dBetaBinom_s} is used when \code{shape1} and \code{shape2} are scalars.
#' In both cases, \code{x} is a vector.
#'
#' @seealso For beta binomial N-mixture models, see \code{\link{dNmixture}}.
#'   For documentation on the beta function, use \code{?stats::dbeta}
#'
#' @examples
#' # Calculate a beta binomial probability with different shape1 and shape2 for each x[i]
#' dBetaBinom_v(x = c(4, 0, 0, 3), N = 10,
#'   shape1 = c(0.5, 0.5, 0.3, 0.5), shape2 = c(0.2, 0.4, 1, 1.2))
#' # or with constant shape1 and shape2
#' dBetaBinom_s(x = c(4, 0, 0, 3), N = 10, shape1 = 0.5, shape2 = 0.5, log = TRUE)


##### Beta binomial support functions #####
#' @rdname dBetaBinom
#' @export
nimBetaFun <- nimbleFunction(
  run = function(a = double(0),
                 b = double(0),
                 log = logical(0)) {
    log_ans <- lgamma(a) + lgamma(b) - lgamma(a + b)
    if (log) return(log_ans)
    else return(exp(log_ans))
    returnType(double(0))
  }, buildDerivs=list(run=list()))

#' @rdname dBetaBinom
#' @export
dBetaBinom_v <- nimbleFunction(
  run = function(x = double(1),
                 N = double(0),
                 shape1 = double(1),
                 shape2 = double(1),
                 len = double(),
                 log = integer(0, default = 0)) {
    logprob <- 0
    lgNp1 <- lgamma(N+1)
    for (i in 1:length(x)) {
      xi <- ADbreak(x[i])
      if(!is.na(xi)){
        logprob <- logprob +
          nimBetaFun(a = x[i] + shape1[i], b = N - x[i] + shape2[i], log = TRUE) -
          nimBetaFun(a = shape1[i], b = shape2[i], log = TRUE) +
          lgNp1 - (lgamma(x[i] + 1) + lgamma(N - x[i] + 1))
      }
    }

    if (log) return(logprob)
    return(exp(logprob))
    returnType(double(0))
  },
  buildDerivs = list(run = list(ignore = c('i','xi')))
)

#' @rdname dBetaBinom
#' @export
dBetaBinom_s <- nimbleFunction(
  run = function(x = double(1),
                 N = double(0),
                 shape1 = double(0),
                 shape2 = double(0),
                 len = double(),
                 log = integer(0, default = 0)) {
    logprob <- 0
    lgNp1 <- lgamma(N+1)
    lbs1s2 <- nimBetaFun(a = shape1, b = shape2, log = TRUE)
    for (i in 1:length(x)) {
      xi <- ADbreak(x[i])
      if(!is.na(xi)){
        logprob <- logprob +
          nimBetaFun(a = x[i] + shape1, b = N - x[i] + shape2, log = TRUE) -
          lbs1s2 +
          lgNp1 - (lgamma(x[i] + 1) + lgamma(N - x[i] + 1))
      }
    }
    if (log) return(logprob)
    return(exp(logprob))
    returnType(double(0))
  },
  buildDerivs = list(run=list(ignore = c('i','xi')))
)

#' @rdname dBetaBinom
#' @export
#' @importFrom stats rbeta
rBetaBinom_v <- nimbleFunction(
  run = function(n = double(0),
                 N = double(0),
                 shape1 = double(1),
                 shape2 = double(1),
                 len = double()) {
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
#' @importFrom stats rbeta
rBetaBinom_s <- nimbleFunction(
  run = function(n = double(0),
                 N = double(0),
                 shape1 = double(0),
                 shape2 = double(0),
                 len = double()) {
    p <- numeric(length=len)
    for (i in 1:len) {
      p[i] <- rbeta(1, shape1, shape2)
    }
    x <- rbinom(len, N, p)
    return(x)
    returnType(double(1))
  })
