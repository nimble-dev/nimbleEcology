#' dNmixture distribution for use in \code{nimble} models
#'
#' \code{dNmixture_s} and \code{dNmixture_v} provide Poisson-Binomial mixture
#' distributions of abundance ("N-mixture") for use in \code{nimble} models.
#' Overdispersion alternatives using the negative binomial distribution (for
#' the abundance submodel) and the beta binomial distribution (for the detection
#' submodel) are also provided.
#'
#' @name dNmixture
#'
#' @aliases dNmixture_s dNmixture_v rNmixture_s rNmixture_v dNmixture_BNB_v
#'   dNmixture_BNB_s dNmixture_BNB_oneObs dNmixture_BBP_v dNmixture_BBP_s
#'   dNmixture_BBP_oneObs dNmixture_BBNB_v dNmixture_BBNB_s
#'   dNmixture_BBNB_oneObs rNmixture_BNB_v rNmixture_BNB_s rNmixture_BNB_oneObs
#'   rNmixture_BBP_v rNmixture_BBP_s rNmixture_BBP_oneObs rNmixture_BBNB_v
#'   rNmixture_BBNB_s rNmixture_BBNB_oneObs
#'
#' @author Ben Goldstein, Lauren Ponisio, and Perry de Valpine
#'
#' @param x vector of integer counts from a series of sampling occasions.
#' @param lambda expected value of the Poisson distribution of true abundance
#' @param theta abundance overdispersion parameter required for negative
#'   binomial (*NB) N-mixture models. The negative binomial is parameterized
#'   such that variance of x is \code{lambda^2 * theta + lambda}
#' @param prob detection probability (scalar for \code{dNmixture_s}, vector for
#'   \code{dNmixture_v}).
#' @param s detection overdispersion parameter required for beta binomial (BB*)
#'   N-mixture models. The beta binomial is parameterized such that variance of
#'   x is \code{V(x) = N * prob * (1-prob) * (N + s) / (s + 1)}
#' @param Nmin minimum abundance to sum over for the mixture probability. Set to
#'   -1 to select automatically (not available for beta binomial variations; see
#'   Details).
#' @param Nmax maximum abundance to sum over for the mixture probability. Set to
#'   -1 to select automatically (not available for beta binomial variations; see
#'   Details).
#' @param len The length of the x vector
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return
#'   probability.
#' @param n number of random draws, each returning a vector of length
#'   \code{len}. Currently only \code{n = 1} is supported, but the argument
#'   exists for standardization of "\code{r}" functions.
#'
#' @details These nimbleFunctions provide distributions that can be
#'     used directly in R or in \code{nimble} hierarchical models (via
#'     \code{\link[nimble]{nimbleCode}} and
#'     \code{\link[nimble]{nimbleModel}}).
#'
#' An N-mixture model defines a distribution for multiple counts (typically of
#' animals, typically made at a sequence of visits to the same site).  The
#' latent number of animals available to be counted, N, follows a Poisson or
#' negative binomial distribution. Each count, \code{x[i]} for visit \code{i},
#' follows a binomial or beta-binomial distribution. The N-mixture distributions
#' calculate the marginal probability of observed counts by summing over the
#' range of latent abundance values.
#'
#' The basic N-mixture model uses Poisson latent abundance with mean
#' \code{lambda} and binomial observed counts  with size (number of trials) N
#' and probability of success (being counted) \code{prob[i]}. This distribution
#' is available in two forms, \code{dNmixture_s} and \code{dNmixture_v}. With
#' \code{dNmixture_s}, detection probability is a scalar, independent of visit,
#' so \code{prob[i]} should be replaced with \code{prob} above.  With
#' \code{dNmixture_v}, detection probability is a vector, with one element for
#' each visit, as written above.
#'
#' We also provide three important variations on the traditional N-mixture
#' model: \code{dNmixture_BNB}, \code{dNmixture_BBP}, and \code{dNmixture_BBNB}.
#' These distributions allow you to replace the Poisson (P) abundance
#' distribution with the negative binomial (NB) and the binomial (B) detection
#' distribution with the beta binomial (BB).
#'
#' Binomial-negative binomial: BNB N-mixture models use a binomial distribution
#' for detection and a negative binomial distribution for abundance with scalar
#' overdispersion parameter \code{theta} (0-Inf). We parameterize such that the
#' variance of the negative binomial is \code{lambda^2 * theta + lambda}, so
#' large \code{theta} indicates a large amount of overdisperison in abundance.
#' The BNB is available in three suffixed forms: \code{dNmixture_BNB_v} is used
#' if \code{prob} varies between observations, \code{dNmixture_BNB_s} is used if
#' \code{prob} is scalar (constant across observations), and
#' \code{dNmixture_BNB_oneObs} is used if only one observation is available at
#' the site (so both x and prob are scalar).
#'
#' Beta-binomial-Poisson: BBP N-mixture uses a beta binomial distribution for
#' detection probabilities and a Poisson distribution for abundance. The beta
#' binomial distribution has scalar overdispersion parameter s (0-Inf). We
#' parameterize such that the variance of the beta binomial is \code{N * prob
#' * (1-prob) * (N + s) / (s + 1)}, with greater s indicating less variance
#' (greater-than-binomial relatedness between observations at the site) and s ->
#' 0 indicating the binomial. The BBP is available in three suffixed forms:
#' \code{dNmixture_BBP_v} is used if \code{prob} varies between observations,
#' \code{dNmixture_BBP_s} is used if \code{prob} is scalar (constant across
#' observations), and \code{dNmixture_BBP_oneObs} is used if only one
#' observation is available at the site (so both x and prob are scalar).
#'
#' Beta-binomial-negative-binomial: dNmixture_BBNB is available using a negative
#' binomial abundance distribution and a beta binomial detection distribution.
#' \code{dNmixture_BBNB} is available with \code{_s}, \code{_v}, and
#' \code{_oneObs} suffixes as above and requires both arguments \code{s} and
#' \code{theta} as parameterized above.
#'
#' The distribution dNmixture_oneObs is not provided as the probability given
#' by the traditional N-mixture distribution for \code{length(x) = 1} is
#' equivalent to \code{dpois(prob * lambda)}.
#'
#' For more explanation, see package vignette
#' (\code{vignette("Introduction_to_nimbleEcology")}).
#'
#' Compared to writing \code{nimble} models with a discrete latent state of
#' abundance N and a separate scalar datum for each count, use of these
#' distributions allows one to directly sum (marginalize) over the discrete
#' latent state N and calculate the probability of all observations for a site
#' jointly.
#'
#' If one knows a reasonable range for summation over possible values of N, the
#' start and end of the range can be provided as \code{Nmin} and \code{Nmax}.
#' Otherwise one can set both to -1, in which case values for \code{Nmin} and
#' \code{Nmax} will be chosen based on the 0.0001 and 0.9999 quantiles of the
#' marginal distributions of each count, using the minimum over counts of the
#' former and the maximum over counts of the latter.
#'
#' The summation over N uses the efficient method given by Meehan et al. (2020,
#' see Appendix B) for the basic Poisson-Binomial case, extended for the
#' overdispersion cases in Goldstein and de Valpine (2022).
#'
#' These are \code{nimbleFunction}s written in the format of user-defined
#' distributions for NIMBLE's extension of the BUGS model language. More
#' information can be found in the NIMBLE User Manual at
#' \href{https://r-nimble.org}{https://r-nimble.org}.
#'
#' When using these distributions in a \code{nimble} model, the left-hand side
#' will be used as \code{x}, and the user should not provide the \code{log}
#' argument.
#'
#' For example, in \code{nimble} model code,
#'
#' \code{observedCounts[i, 1:T] ~ dNmixture_v(lambda[i],
#' prob[i, 1:T],
#' Nmin, Nmax, T)}
#'
#' declares that the \code{observedCounts[i, 1:T]} (observed counts
#' for site \code{i}, for example) vector follows an N-mixture
#' distribution with parameters as indicated, assuming all the
#' parameters have been declared elsewhere in the model. As above,
#' \code{lambda[i]} is the mean of the abundance distribution at site
#' i, \code{prob[i, 1:T]} is a vector of detection probabilities at
#' site i, and \code{T} is the number of observation occasions. This
#' will invoke (something like) the following call to
#' \code{dNmixture_v} when \code{nimble} uses the model such as for
#' MCMC:
#'
#' \code{dNmixture_v(observedCounts[i, 1:T], lambda[i],
#' prob[i, 1:T],
#' Nmin, Nmax, T, log = TRUE)}
#'
#' If an algorithm using a \code{nimble} model with this declaration
#' needs to generate a random draw for \code{observedCounts[1:T]}, it
#' will make a similar invocation of \code{rNmixture_v}, with \code{n = 1}.
#'
#' If the observation probabilities are visit-independent, one would use:
#'
#' \code{observedCounts[1:T] ~ dNmixture_s(observedCounts[i, 1:T], lambda[i],
#' prob[i],
#' Nmin, Nmax, T)}
#'
#' @section Notes for use with automatic differentiation:
#'
#' The N-mixture distributions are the only ones in \code{nimbleEcology} for which
#' one must use different versions when AD support is needed. See
#' \code{\link{dNmixtureAD}}.
#'
#' @return
#' For \code{dNmixture_s} and \code{dNmixture_v}: the probability (or likelihood) or log
#' probability of observation vector \code{x}.
#'
#' For \code{rNmixture_s} and \code{rNmixture_v}: a simulated detection history, \code{x}.
#'
#' @seealso For occupancy models dealing with detection/nondetection data,
#' see \code{\link{dOcc}}.
#' For dynamic occupancy, see \code{\link{dDynOcc}}.
#'
#' @references
#'
#' D. Turek, P. de Valpine and C. J. Paciorek. 2016. Efficient Markov chain Monte
#' Carlo sampling for hierarchical hidden Markov models. Environmental and
#' Ecological Statistics 23:549–564. DOI 10.1007/s10651-016-0353-z
#'
#' Meehan, T. D., Michel, N. L., & Rue, H. 2020. Estimating Animal Abundance
#' with N-Mixture Models Using the R—INLA Package for R. Journal of Statistical
#' Software, 95(2). https://doi.org/10.18637/jss.v095.i02
#'
#' Goldstein, B.R., and P. de Valpine. 2022. Comparing N-mixture Models and
#' GLMMs for Relative Abundance Estimation in a Citizen Science Dataset.
#' Scientific Reports 12: 12276. DOI:10.1038/s41598-022-16368-z
#'
#' @examples
#' # Set up constants and initial values for defining the model
#' len <- 5 # length of dataset
#' dat <- c(1,2,0,1,5) # A vector of observations
#' lambda <- 10 # mean abundance
#' prob <- c(0.2, 0.3, 0.2, 0.1, 0.4) # A vector of detection probabilities
#'
#' # Define code for a nimbleModel
#'  nc <- nimbleCode({
#'    x[1:5] ~ dNmixture_v(lambda, prob = prob[1:5],
#'                         Nmin = -1, Nmax = -1, len = 5)
#'
#'    lambda ~ dunif(0, 1000)
#'
#'    for (i in 1:5) {
#'      prob[i] ~ dunif(0, 1)
#'    }
#'  })
#'
#' # Build the model
#' nmix <- nimbleModel(nc,
#'                     data = list(x = dat),
#'                     inits = list(lambda = lambda,
#'                                  prob = prob))
#' # Calculate log probability of data from the model
#' nmix$calculate()
#' # Use the model for a variety of other purposes...

# nimbleOptions(checkNimbleFunction = FALSE)


##### Regular N-mixture #####
NULL
#' @rdname dNmixture
#' @importFrom stats qpois
#' @export
dNmixture_v <- nimbleFunction(
    run = function(x = double(1),
                   lambda = double(),
                   prob = double(1),
                   Nmin = integer(0, default = -1),
                   Nmax = integer(0, default = -1),
                   len = integer(),
                   log = integer(0, default = 0)) {
  if (length(x) != len) stop("in dNmixture_v, len must equal length(x).")
  if (len != length(prob)) stop("in dNmixture_v, len must equal length(prob).")
  if (lambda < 0)
    if (log) return(-Inf) else return(0)
  ## For each x, the conditional distribution of (N - x | x) is pois(lambda * (1-p))
  ## We determine the lowest N and highest N at extreme quantiles and sum over those.
  if (Nmin == -1)
    Nmin <- min(x + qpois(0.00001, lambda * (1 - prob)))
  if (Nmax == -1)
    Nmax <- max(x + qpois(0.99999, lambda * (1 - prob)))

  max_x <- 0
  for (i in 1:length(x)){
    if(!is.na(x[i])){
      if(x[i] > max_x) max_x <- x[i]
    }
  }
  Nmin <- max( max_x, Nmin ) ## set Nmin to at least the largest x

  sum_log_dbinom <- 0
  sum_log_one_m_prob <- 0
  any_not_na <- FALSE
  for (i in 1:length(x)){
    if(!is.na(x[i])){
      sum_log_one_m_prob <- sum_log_one_m_prob + log(1 - prob[i])
      sum_log_dbinom <- sum_log_dbinom + dbinom(x[i], size = Nmin, prob = prob[i], log=TRUE)
      any_not_na <- TRUE
    }
  }

  logProb <- 0
  if(any_not_na){
    logProb <- dNmixture_steps(x, lambda, Nmin, Nmax, sum_log_one_m_prob,
                               sum_log_dbinom)
  }
  if (log) return(logProb)
  else return(exp(logProb))
  returnType(double())
}
)

NULL
#' @rdname dNmixture
#' @importFrom stats qpois
#' @export
dNmixture_s <- nimbleFunction(
    run = function(x = double(1),
                   lambda = double(),
                   prob = double(),
                   Nmin = double(0, default = -1),
                   Nmax = double(0, default = -1),
                   len = double(),
                   log = integer(0, default = 0)) {
  if (length(x) != len) stop("in dNmixture_s, len must equal length(x).")

  # Lambda cannot be negative
  if (lambda < 0)
    if (log) return(-Inf) else return(0)
  ## For each x, the conditional distribution of (N - x | x) is pois(lambda * (1-p))
  ## We determine the lowest N and highest N at extreme quantiles and sum over those.
  if (Nmin == -1)
    Nmin <- min(x + qpois(0.00001, lambda * (1 - prob)))
  if (Nmax == -1)
    Nmax <- max(x + qpois(0.99999, lambda * (1 - prob)))
  
  max_x <- 0
  for (i in 1:length(x)){
    if(!is.na(x[i])){
      if(x[i] > max_x) max_x <- x[i]
    }
  }
  Nmin <- max( max_x, Nmin ) ## set Nmin to at least the largest x

  sum_log_dbinom <- 0
  sum_log_one_m_prob <- 0
  any_not_na <- FALSE
  for (i in 1:length(x)){
    if(!is.na(x[i])){
      sum_log_one_m_prob <- sum_log_one_m_prob + log(1 - prob)
      sum_log_dbinom <- sum_log_dbinom + dbinom(x[i], size = Nmin, prob = prob, log=TRUE)
      any_not_na <- TRUE
    }
  }

  logProb <- 0
  if(any_not_na){
    logProb <- dNmixture_steps(x, lambda, Nmin, Nmax, sum_log_one_m_prob,
                               sum_log_dbinom)
  }
  if (log) return(logProb)
  else return(exp(logProb))
  returnType(double())
 }
)

NULL
#' @rdname dNmixture
#' @importFrom stats rpois rbinom
#' @export
rNmixture_v <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 prob = double(1),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double()) {
    if (n != 1) stop("rNmixture_v only works for n = 1")
    if (length(prob) != len) stop("In rNmixture_v, len must equal length(prob).")
    trueN <- rpois(1, lambda)
    ans <- numeric(len)
    for (i in 1:len)
      ans[i] <- rbinom(n = 1, size = trueN, prob = prob[i])
    return(ans)
    returnType(double(1))
})

NULL
#' @rdname dNmixture
#' @importFrom stats rpois rbinom
#' @export
rNmixture_s <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 prob = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double()) {
    if (n != 1) stop("rNmixture_v only works for n = 1")
    trueN <- rpois(1, lambda)
    ans <- numeric(len)
    for (i in 1:len)
      ans[i] <- rbinom(n = 1, size = trueN, prob = prob)
    return(ans)
    returnType(double(1))
})


#' @rdname dNmixture
#' @importFrom stats qnbinom
#' @export
dNmixture_BNB_v <- nimbleFunction(
  run = function(x = double(1),
                 lambda = double(),
                 theta = double(),
                 prob = double(1),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double(),
                 log = integer(0, default = 0)) {
    if (length(x) != len) stop("in dNmixture_BNB_v, len must equal length(x).")
    if (len != length(prob)) stop("in dNmixture_BNB_v, len must equal length(prob).")

    if (theta <= 0)
      if (log) return(-Inf) else return(0)
    if (lambda < 0)
      if (log) return(-Inf) else return(0)
    ## For each x, the conditional distribution of (N - x | x) is
    ## a negative binomial with overdispersion parameter (theta / (1 + y * theta))
    ## and mean omega / ((theta / (1 + y * theta)) * (1 - omega)) where
    ## omega = (1 - p) * (theta * lambda / (1 + theta * lambda))
    ## We determine the lowest N and highest N at extreme quantiles and sum over those.
    theta_cond <- theta / (1 + x * theta)
    omega <- (1 - prob) * (theta * lambda / (1 + theta * lambda))
    lambda_cond <- omega / (theta_cond * (1 - omega))
    r_cond <- 1 / theta_cond
    pNB_cond <- 1 / (1 + theta_cond * lambda_cond)

    if (Nmax == -1){
      Nmax <- 0
      for (i in 1:length(x)){
        if(!is.na(x[i])){
          Nmax_cand <- x[i] + qnbinom(0.99999, size = r_cond[i], prob = pNB_cond[i])
          if(Nmax_cand > Nmax) Nmax <- Nmax_cand
        }
      }
    }

    if(Nmin == -1){
      Nmin <- Nmax
      for (i in 1:length(x)){
        if(!is.na(x[i])){
          Nmin_cand <- x[i] + qnbinom(0.00001, size = r_cond[i], prob = pNB_cond[i])
          if(Nmin_cand < Nmin) Nmin <- Nmin_cand
        }
      }
    }

    max_x <- 0
    for (i in 1:length(x)){
      if(!is.na(x[i])){
        if(x[i] > max_x) max_x <- x[i]
      }
    }

    Nmin <- max( max_x, Nmin ) ## set Nmin to at least the largest x

    sum_log_dbinom <- 0
    sum_log_one_m_prob <- 0
    any_not_na <- FALSE
    for (i in 1:length(x)){
      if(!is.na(x[i])){
        sum_log_one_m_prob <- sum_log_one_m_prob + log(1 - prob[i])
        sum_log_dbinom <- sum_log_dbinom + dbinom(x[i], size = Nmin, prob = prob[i], log=TRUE)
        any_not_na <- TRUE
      }
    }

    logProb <- 0
    if(any_not_na){
      logProb <- dNmixture_BNB_steps(x,lambda,theta,Nmin,Nmax, sum_log_one_m_prob,
                                     sum_log_dbinom)
    }
    if (log) return(logProb)
    else return(exp(logProb))
    returnType(double())
  }
)

##### dNmixture_BNB_s #####
NULL
#' @rdname dNmixture
#' @importFrom stats qnbinom
#' @export
dNmixture_BNB_s <- nimbleFunction(
  run = function(x = double(1),
                 lambda = double(),
                 theta = double(),
                 prob = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double(),
                 log = integer(0, default = 0)) {
    if (length(x) != len) stop("in dNmixture_BNB_s, len must equal length(x).")
    if (theta <= 0)
      if (log) return(-Inf) else return(0)
    if (lambda < 0)
      if (log) return(-Inf) else return(0)
    ## See above for comments
    theta_cond <- theta / (1 + x * theta)
    omega <- (1 - prob) * (theta * lambda / (1 + theta * lambda))
    lambda_cond <- omega / (theta_cond * (1 - omega))
    r_cond <- 1 / theta_cond
    pNB_cond <- 1 / (1 + theta_cond * lambda_cond)

    if (Nmax == -1){
      Nmax <- 0
      for (i in 1:length(x)){
        if(!is.na(x[i])){
          Nmax_cand <- x[i] + qnbinom(0.99999, size = r_cond[i], prob = pNB_cond[i])
          if(Nmax_cand > Nmax) Nmax <- Nmax_cand
        }
      }
    }

    if(Nmin == -1){
      Nmin <- Nmax
      for (i in 1:length(x)){
        if(!is.na(x[i])){
          Nmin_cand <- x[i] + qnbinom(0.00001, size = r_cond[i], prob = pNB_cond[i])
          if(Nmin_cand < Nmin) Nmin <- Nmin_cand
        }
      }
    }

    max_x <- 0
    for (i in 1:length(x)){
      if(!is.na(x[i])){
        if(x[i] > max_x) max_x <- x[i]
      }
    }

    Nmin <- max( max_x, Nmin ) ## set Nmin to at least the largest x

    sum_log_dbinom <- 0
    sum_log_one_m_prob <- 0
    any_not_na <- FALSE
    for (i in 1:length(x)){
      if(!is.na(x[i])){
        sum_log_one_m_prob <- sum_log_one_m_prob + log(1 - prob)
        sum_log_dbinom <- sum_log_dbinom + dbinom(x[i], size = Nmin, prob = prob, log=TRUE)
        any_not_na <- TRUE
      }
    }

    logProb <- 0
    if(any_not_na){
      logProb <- dNmixture_BNB_steps(x,lambda,theta,Nmin,Nmax, sum_log_one_m_prob,
                                     sum_log_dbinom)
    }
    if (log) return(logProb)
    else return(exp(logProb))
    returnType(double())
  }
)

##### dNmixture_BNB_oneObs #####
NULL
#' @rdname dNmixture
#' @export
dNmixture_BNB_oneObs <- nimbleFunction(
  run = function(x = double(),
                 lambda = double(),
                 theta = double(),
                 prob = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 log = integer(0, default = 0)) {
    xvec <- numeric(value = x, length=1)
    return(dNmixture_BNB_s(xvec,lambda,theta,prob,Nmin,Nmax,1,log))
    returnType(double())
  }
)

#' @rdname dNmixture
#' @export
dNmixture_BBP_v <- nimbleFunction(
  run = function(x = double(1),
                 lambda = double(),
                 prob = double(1),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double(),
                 log = integer(0, default = 0)) {
    if (length(x) != len) stop("in dNmixture_BBP_v, len must equal length(x).")
    if (len != length(prob)) stop("in dNmixture_BBP_v, len must equal length(prob).")
    if (s <= 0)
      if (log) return(-Inf) else return(0)
    alpha <- prob * s
    beta <- s - prob * s
    if (lambda < 0)
      if (log) return(-Inf) else return(0)
    ## For beta-binomial N-mixtures , the conditional distribution of (N - x |
    ## x) doesn't have a nice closed-form expression.
    if (Nmin == -1 | Nmax == -1) {
      stop("Dynamic choice of Nmin/Nmax is not supported for beta binomial N-mixtures.")
    }
    Nmin <- max( max(x), Nmin ) ## set Nmin to at least the largest x
    logProb <- dNmixture_BBP_steps(x, beta-x, lambda, s, Nmin, Nmax,
                                   dBetaBinom_v(x, Nmin, alpha, beta, len, log = TRUE))
    if (log) return(logProb)
    else return(exp(logProb))
    returnType(double())
  }
)

##### dNmixture_BBP_s #####
NULL
#' @rdname dNmixture
#' @export
dNmixture_BBP_s <- nimbleFunction(
  run = function(x = double(1),
                 lambda = double(),
                 prob = double(),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double(),
                 log = integer(0, default = 0)) {
    if (length(x) != len) stop("in dNmixture_BBP_s, len must equal length(x).")
    if (s <= 0)
      if (log) return(-Inf) else return(0)
    alpha <- prob * s
    beta <- s - prob * s
    if (lambda < 0)
      if (log) return(-Inf) else return(0)
    ## For beta-binomial N-mixtures , the conditional distribution of (N - x |
    ## x) doesn't have a nice closed-form expression.
    if (Nmin == -1 | Nmax == -1) {
      stop("Dynamic choice of Nmin/Nmax is not supported for beta binomial N-mixtures.")
    }
    Nmin <- max( max(x), Nmin ) ## set Nmin to at least the largest x
    logProb <- dNmixture_BBP_steps(x, beta-x, lambda, s, Nmin, Nmax,
                                   dBetaBinom_s(x, Nmin, alpha, beta, len, log = TRUE))
    if (log) return(logProb)
    else return(exp(logProb))
    returnType(double())
  }
)

##### dNmixture_BBP_oneObs #####
NULL
#' @rdname dNmixture
#' @export
dNmixture_BBP_oneObs <- nimbleFunction(
  run = function(x = double(),
                 lambda = double(),
                 prob = double(),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 log = integer(0, default = 0)) {
    xvec <- numeric(value = x, length = 1)
    return(dNmixture_BBP_s(xvec,lambda,prob,s,Nmin,Nmax,1,log))
    returnType(double())
  }
)

#' @rdname dNmixture
#' @export
dNmixture_BBNB_v <- nimbleFunction(
  run = function(x = double(1),
                 lambda = double(),
                 theta = double(),
                 prob = double(1),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double(),
                 log = integer(0, default = 0)) {
    if (length(x) != len) stop("in dNmixture_BBNB_v, len must equal length(x).")
    if (len != length(prob)) stop("in dNmixture_BBNB_v, len must equal length(prob).")
    if (s <= 0)
      if (log) return(-Inf) else return(0)
    if (theta <= 0)
      if (log) return(-Inf) else return(0)
    alpha <- prob * s
    beta <- s - prob * s
    if (lambda < 0)
      if (log) return(-Inf) else return(0)
    ## see comments above
    if (Nmin == -1 | Nmax == -1) {
      stop("Dynamic choice of Nmin/Nmax is not supported for beta binomial N-mixtures.")
    }
    Nmin <- max( max(x), Nmin ) ## set Nmin to at least the largest x
    logProb <- dNmixture_BBNB_steps(x, beta-x,lambda,theta,s,Nmin,Nmax,
                                    dBetaBinom_v(x, Nmin, alpha, beta, len, log = TRUE))
    if (log) return(logProb)
    else return(exp(logProb))
    returnType(double())
  }
)

##### dNmixture_BBNB_s #####
NULL
#' @rdname dNmixture
#' @export
dNmixture_BBNB_s <- nimbleFunction(
  run = function(x = double(1),
                 lambda = double(),
                 theta = double(),
                 prob = double(),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double(),
                 log = integer(0, default = 0)) {
    if (length(x) != len) stop("in dNmixture_BBNB_s, len must equal length(x).")
    if (s <= 0)
      if (log) return(-Inf) else return(0)
    if (theta <= 0)
      if (log) return(-Inf) else return(0)
    r <- 1 / theta
    pNB <- 1 / (1 + theta * lambda)
    alpha <- prob * s
    beta <- s - prob * s
    if (lambda < 0)
      if (log) return(-Inf) else return(0)
    ## See comments above
    if (Nmin == -1 | Nmax == -1) {
      stop("Dynamic choice of Nmin/Nmax is not supported for beta binomial N-mixtures.")
    }
    Nmin <- max( max(x), Nmin ) ## set Nmin to at least the largest x
    logProb <- dNmixture_BBNB_steps(x, beta-x,lambda,theta,s,Nmin,Nmax,
                                    dBetaBinom_s(x, Nmin, alpha, beta, len, log = TRUE))
    if (log) return(logProb)
    else return(exp(logProb))
    returnType(double())
  }
)

##### dNmixture_BBNB_oneObs #####
NULL
#' @rdname dNmixture
#' @export
dNmixture_BBNB_oneObs <- nimbleFunction(
  run = function(x = double(),
                 lambda = double(),
                 theta = double(),
                 prob = double(),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 log = integer(0, default = 0)) {
    xvec <- numeric(x, length = 1)
    return(dNmixture_BBNB_s(xvec, lambda, theta, prob, s, Nmin, Nmax, 1, log))
    returnType(double())
  }
)

##### rNmixture extensions #####
NULL
#' @rdname dNmixture
#' @importFrom stats rnbinom
#' @export
rNmixture_BNB_v <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 theta = double(),
                 prob = double(1),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double()) {

    if (n != 1) stop("rNmixture* only works for n = 1")
    if (length(prob) != len) stop("In rNmixture*, len must equal length(prob).")

    r <- 1 / theta
    p <- 1 / (1 + theta * lambda)

    trueN <- rnbinom(1, size = r, prob = p)
    ans <- numeric(len)
    for (i in 1:len) {
      ans[i] <- rbinom(n = 1, size = trueN, prob = prob[i])
    }

    return(ans)
    returnType(double(1))
  })

NULL
#' @rdname dNmixture
#' @importFrom stats rnbinom
#' @export
rNmixture_BNB_s <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 theta = double(),
                 prob = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double()) {

    if (n != 1) stop("rNmixture* only works for n = 1")

    r <- 1 / theta
    p <- 1 / (1 + theta * lambda)

    trueN <- rnbinom(1, size = r, prob = p)
    ans <- numeric(len)
    for (i in 1:len) {
      ans[i] <- rbinom(n = 1, size = trueN, prob = prob)
    }

    return(ans)
    returnType(double(1))
  })


NULL
#' @rdname dNmixture
#' @importFrom stats rnbinom
#' @export
rNmixture_BNB_oneObs <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 theta = double(),
                 prob = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1)) {
    if (n != 1) stop("rNmixture* only works for n = 1")

    r <- 1 / theta
    p <- 1 / (1 + theta * lambda)
    trueN <- rnbinom(1, size = r, prob = p)

    ans <- rbinom(n = 1, size = trueN, prob = prob)

    return(ans)
    returnType(double())
  })


NULL
#' @rdname dNmixture
#' @importFrom stats rpois
#' @export
rNmixture_BBP_v <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 prob = double(1),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double()) {
    if (n != 1) stop("rNmixture* only works for n = 1")
    if (length(prob) != len) stop("In rNmixture*, len must equal length(prob).")

    alpha <- prob * s
    beta <- s - prob * s

    trueN <- rpois(1, lambda = lambda)
    ans <- rBetaBinom_v(n = 1, N = trueN, shape1 = alpha, shape2 = beta, len = len)

    return(ans)
    returnType(double(1))
  })

NULL
#' @rdname dNmixture
#' @importFrom stats rpois
#' @export
rNmixture_BBP_s <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 prob = double(),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double()) {
    if (n != 1) stop("rNmixture* only works for n = 1")

    alpha <- prob * s
    beta <- s - prob * s

    trueN <- rpois(1, lambda = lambda)
    ans <- rBetaBinom_s(n = 1, N = trueN,
                      shape1 = alpha, shape2 = beta, len = len)

    return(ans)
    returnType(double(1))
  })

NULL
#' @rdname dNmixture
#' @importFrom stats rpois
#' @export
rNmixture_BBP_oneObs <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 prob = double(),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1)) {
    if (n != 1) stop("rNmixture* only works for n = 1")

    alpha <- prob * s
    beta <- s - prob * s

    trueN <- rpois(1, lambda = lambda)
    ans <- rBetaBinom_s(n = 1, N = trueN, shape1 = alpha, shape2 = beta, len = 1)

    return(ans[1])
    returnType(double())
  })

NULL
#' @rdname dNmixture
#' @importFrom stats rnbinom
#' @export
rNmixture_BBNB_v <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 theta = double(),
                 prob = double(1),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double()) {
    if (n != 1) stop("rNmixture* only works for n = 1")

    alpha <- prob * s
    beta <- s - prob * s
    r <- 1 / theta
    p <- 1 / (1 + theta * lambda)

    trueN <- rnbinom(1, size = r, prob = p)
    ans <- rBetaBinom_v(n = 1, N = trueN, shape1 = alpha, shape2 = beta, len = len)

    return(ans)
    returnType(double(1))
  })

NULL
#' @rdname dNmixture
#' @importFrom stats rnbinom
#' @export
rNmixture_BBNB_s <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 theta = double(),
                 prob = double(),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double()) {
    if (n != 1) stop("rNmixture* only works for n = 1")

    alpha <- prob * s
    beta <- s - prob * s
    r <- 1 / theta
    p <- 1 / (1 + theta * lambda)

    trueN <- rnbinom(1, size = r, prob = p)
    ans <- rBetaBinom_s(n = 1, N = trueN,
                      shape1 = alpha, shape2 = beta, len = len)

    return(ans)
    returnType(double(1))
  })

NULL
#' @rdname dNmixture
#' @importFrom stats rnbinom
#' @export
rNmixture_BBNB_oneObs <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 theta = double(),
                 prob = double(),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1)) {
    if (n != 1) stop("rNmixture* only works for n = 1")

    alpha <- prob * s
    beta <- s - prob * s
    r <- 1 / theta
    p <- 1 / (1 + theta * lambda)

    trueN <- rnbinom(1, size = r, prob = p)
    ans <- rBetaBinom_s(n = 1, N = trueN, shape1 = alpha, shape2 = beta, len = 1)
    return(ans[1])
    returnType(double())
  })
