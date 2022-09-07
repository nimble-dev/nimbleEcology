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
#' @author Juniper Simonis
#'
#' @param x vector of integer counts from a series of sampling occasions.
#' @param mu expected value of the (negative binomial or Poisson) distribution of true abundance.
#' @param r shape parameter defining overdispersion. As \code{r} approaches 0, the negative binomial converges to a Poisson.
#' @param p detection probability (scalar for \code{dNmixture_M\*_s}, vector for \code{dNmixture_M\*_v}).
#' @param J integer number of searches.
#' @param log \code{TRUE} or \code{1} to return log probability. \code{FALSE} or \code{0} to return probability.
#' @param n number of random draws, each returning a vector of length \code{len}. Currently only \code{n = 1} is supported, but the argument exists for standardization of "\code{r}" functions.
#'
#' @details These nimbleFunctions provide distributions that can be used directly in R or in \code{nimble} hierarchical models (via \code{\link[nimble]{nimbleCode}} and \code{\link[nimble]{nimbleModel}}). \cr
#'          The distributions are implemented in closed-form following Haines (2020), which avoids infinite sum-based calculations.
#'
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
#' #  mu: mean abundance
#' #  r:  scale parameter on abundance distribution
#' #  p:  search-specific detection probabilities
#'
#'  J <- 10
#'  mu <- 5
#'  r  <- 2
#'  p  <- runif(J, 0.4, 0.7)
#'
#'  mut <- log(mu)
#'  rt  <- log(r)
#'
#' # Generate data
#'
#'   yv <- rNmixture_MNB_v(n = 1, mu = mu, p = p, r = r, J = J)
#'
#' # Write the model code
#'
#'   nc <- nimbleCode({
#'
#'     mut ~  dnorm(0, 1/2)
#'     mu  <- exp(mut)
#'
#'     rt ~  dnorm(0, 0.1)
#'     r  <- exp(rt)
#'
#'     for (j in 1:J) {
#'       p[j] ~ dunif(0, 1)
#'     }
#'
#'     x[1:J] ~ dNmixture_MNB_v(mu = mu, p = p[1:J], r = r, J = J)
#'
#'   })
#'
#' # Construct the model object
#'
#'   nmix <- nimbleModel(nc,
#'                       constants = list(J = J),
#'                       data = list(x = yv),
#'                       inits = list(mut = mut,
#'                                    rt  = rt,
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
    run = function(x   = double(1),
                   mu  = double(),
                   p   = double(),
                   r   = double(),
                   J   = double(),
                   log = integer(0, default = 0)) {

  x_tot <- sum(x)
  x_miss <- sum(x * seq(0, J - 1))

  term1   <- lgamma(r + x_tot) - lgamma(r) - sum(lfactorial(x))
  term2   <- r * log(r) + x_tot * log(mu)
  term3   <- x_tot * log(p) + x_miss * log(1 - p)
  term4   <- -(x_tot + r) * log(r + mu * (1 - (1 - p) ^ J))
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
  run = function(n  = integer(),
                 mu = double(),
                 p  = double(),
                 r  = double(),
                 J  = double()) {


    prob <- numeric(J + 1)
    for (i in 1:(J)) {
      prob[i] <- pow(1 - p, i - 1) * p
    }
    prob[J + 1] <- 1 - sum(prob[1:J])

    ans <- numeric(J + 1)
    n <- rnbinom(n = 1, size = r, prob = 1/(1 + (1/r) * mu))
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
    run = function(x   = double(1),
                   mu  = double(),
                   p   = double(1),
                   r   = double(),
                   J   = double(),
                   log = integer(0, default = 0)) {

  x_tot <- sum(x)
  x_miss <- sum(x * seq(0, J - 1))
  pp <- c(0, p)
  prob <- numeric(J)
  for (j in 1:J) {
    prob[j] <- prod(1 - pp[1:j]) * p[j]
  }
  ptot <- sum(prob)

  # from _s:
  # term1   <- lgamma(r + x_tot) - lgamma(r) - sum(lfactorial(x))
  # term2   <- r * log(r) + x_tot * log(mu)
  # term3   <- x_tot * log(p) + x_miss * log(1 - p)
  # term4   <- -(x_tot + r) * log(r + mu * (1 - (1 - p) ^ J))
  #

  term1   <- lgamma(r + x_tot) - lgamma(r) - sum(lfactorial(x))
  term2   <- r * log(r) + x_tot * log(mu)
  term3   <- sum(x * log(prob))
  term4   <- -(x_tot + r) * log(r + mu * ptot)
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
  run = function(n  = integer(),
                 mu = double(),
                 p  = double(1),
                 r  = double(),
                 J  = double()) {


    pp <- c(0, p)
    prob <- numeric(J + 1)
    for (j in 1:J) {
      prob[j] <- prod(1 - pp[1:j]) * p[j]
    }
    prob[J + 1] <- 1 - sum(prob[1:J])

    ans <- numeric(J + 1)
    n <- rnbinom(n = 1, size = r, prob = 1/(1 + (1/r) * mu))
    if (n > 0) {
      ans <- rmulti(n = 1, size = n, prob = prob)
    }

    return(ans[1:J])
    returnType(double(1))
})
