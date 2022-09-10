#
# Multinomial Negative Binomial mixture models with scalar and vector p
#


# dNmixture_M
#' @title Multinomial N-mixture distribution for use in \code{nimble} models
#'
#' @description \code{dNmixture_MNB_s} and \code{dNmixture_MNB_v} provide Multinomial Negative Binomial (MNB) mixture distributions of abundance ("N-mixture") for use in \code{nimble} models.
#'
#' @name dNmixture_M
#' @aliases dNmixture_MNB_s dNmixture_MNB_v rNmixture_MNB_s rNmixture_MNB_v
#' dNmixture_MP_s dNmixture_MP_v rNmixture_MP_s rNmixture_MP_v
#'
#' @author Juniper Simonis, Ben Goldstein
#'
#' @param x vector of integer counts from a series of sampling occasions.
#' @param lambda expected value of the (negative binomial or Poisson)
#'        distribution of true abundance.
#' @param theta overdispersion parameter required for negative binomial
#'        (*NB) models. theta is parameterized such that variance of
#'        the negative binomial variable x is \code{lambda^2 * theta + lambda}
#' @param p detection probability (scalar for \code{dNmixture_M\*_s},
#'        vector for \code{dNmixture_M\*_v}). If a vector, p[i] is the probability
#'        that an individual is observed in observation category i and no other
#'        categories. \code{p} must sum to less than or equal to 1, although this
#'        is not checked by the function for computational efficiency. The
#'        probability of an individual going unobserved is \code{1 - sum(p)}.
#' @param J integer number of searches.
#' @param log \code{TRUE} or \code{1} to return log probability. \code{FALSE} or
#'        \code{0} to return probability.
#' @param n number of random draws, each returning a vector of length
#'        \code{len}. Currently only \code{n = 1} is supported, but the argument
#'        exists for standardization of "\code{r}" functions.
#'
#' @details These nimbleFunctions provide distributions that can be used directly in R or in \code{nimble} hierarchical models (via \code{\link[nimble]{nimbleCode}} and \code{\link[nimble]{nimbleModel}}). \cr
#'          The distributions are implemented in closed-form following Haines (2020), which avoids infinite sum-based calculations.
#'
#' @return For \code{dNmixture_s} and \code{dNmixture_v}: the probability (or likelihood) or log probability of observation vector \code{x}. \cr
#'         For \code{rNmixture_s} and \code{rNmixture_v}: a simulated detection history, \code{x}.
#'
#' @seealso For binomial N-mixture models, see \code{\link{dNmixture}}. \cr
#'          For occupancy models dealing with detection/nondetection data, see \code{\link{dOcc}}. \cr
#'          For dynamic occupancy, see \code{\link{dDynOcc}}.
#'
#' @references
#'   L. M. Haines. 2020. Multinomial N-mixture models for removal sampling. Biometrics 76:540-548. DOI 10.1111/biom.13147
#'
#' @examples
#'
#' # Set up variables and parameters
#' #  J:  number of visits
#' #  lambda: mean abundance
#' #  theta:  scale parameter on abundance distribution
#' #  p:  search-specific detection probabilities
#'
#'  J <- 10
#'  lambda <- 5
#'  theta  <- 2
#'  p  <- runif(J, 0.4, 0.7)
#'
#'  lambdat <- log(mu)
#'  thetat  <- log(theta)
#'
#' # Generate data
#'
#'   yv <- rNmixture_MNB_v(n = 1, lambda = lambda, p = p, theta = theta, J = J)
#'
#' # Write the model code
#'
#'   nc <- nimbleCode({
#'
#'     lambdat ~  dnorm(0, 1/2)
#'     lambda  <- exp(lambdat)
#'
#'     thetat ~  dnorm(0, 0.1)
#'     theta  <- exp(thetat)
#'
#'     for (j in 1:J) {
#'       p[j] ~ dunif(0, 1)
#'     }
#'
#'     x[1:J] ~ dNmixture_MNB_v(lambda = lambda, p = p[1:J], theta = theta, J = J)
#'
#'   })
#'
#' # Construct the model object
#'
#'   nmix <- nimbleModel(nc,
#'                       constants = list(J = J),
#'                       data = list(x = yv),
#'                       inits = list(lambdat = lambdat,
#'                                    thetat  = thetat,
#'                                    p   = p))
#'
#' # Calculate log probability of data from the model
#'
#'   nmix$calculate()
#'
#'
# nimbleOptions(checkNimbleFunction = FALSE)
NULL

#' @rdname dNmixture_M
#'
#' @export
#'
dNmixture_MNB_s <- nimbleFunction(
    run = function(x       = double(1),
                   lambda  = double(),
                   p       = double(),
                   theta   = double(),
                   J       = double(),
                   log     = integer(0, default = 0)) {

    p_vec <- rep(p, J)

    r <- 1/theta

    x_tot <- sum(x)
    x_miss <- sum(x * seq(0, J - 1))

    ptot <- sum(p_vec)

    term1   <- lgamma(r + x_tot) - lgamma(r) - sum(lfactorial(x))
    term2   <- r * log(r) + x_tot * log(lambda)
    term3   <- sum(x * log(p_vec))
    term4   <- -(x_tot + r) * log(r + lambda * ptot)
    logProb <- term1 + term2 + term3 + term4


  if (log) return(logProb)
  else return(exp(logProb))
  returnType(double())

})

#' @rdname dNmixture_M
#'
#' @export
#'
rNmixture_MNB_s <- nimbleFunction(
  run = function(n      = integer(),
                 lambda = double(),
                 p      = double(),
                 theta  = double(),
                 J      = double()) {


    prob <- c(rep(p, J), 1 - (p*J))

    ans <- numeric(J + 1)
    n <- rnbinom(n = 1, size = 1/theta, prob = 1/(1 + theta * lambda))
    if (n > 0) {
      ans <- rmulti(n = 1, size = n, prob = prob)
    }

    return(ans[1:J])
    returnType(double(1))
})


#' @rdname dNmixture_M
#'
#' @export
#'
dNmixture_MNB_v <- nimbleFunction(
    run = function(x       = double(1),
                   lambda  = double(),
                   p       = double(1),
                   theta   = double(),
                   J       = double(),
                   log     = integer(0, default = 0)) {

  r <- 1/theta

  x_tot <- sum(x)
  x_miss <- sum(x * seq(0, J - 1))

  ptot <- sum(p)

  term1   <- lgamma(r + x_tot) - lgamma(r) - sum(lfactorial(x))
  term2   <- r * log(r) + x_tot * log(lambda)
  term3   <- sum(x * log(p))
  term4   <- -(x_tot + r) * log(r + lambda * ptot)
  logProb <- term1 + term2 + term3 + term4

  if (log) return(logProb)
  else return(exp(logProb))
  returnType(double())

})


#' @rdname dNmixture_M
#'
#' @export
#'
rNmixture_MNB_v <- nimbleFunction(
  run = function(n      = integer(),
                 lambda = double(),
                 p      = double(1),
                 theta  = double(),
                 J      = double()) {

    prob <- c(p, 1 - sum(p))

    ans <- numeric(J + 1)
    n <- rnbinom(n = 1, size = 1/theta, prob = 1/(1 + theta * lambda))
    if (n > 0) {
      ans <- rmulti(n = 1, size = n, prob = prob)
    }

    return(ans[1:J])
    returnType(double(1))
})
