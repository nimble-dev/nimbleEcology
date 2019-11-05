# dNmixture
#' N-mixture distribution for use in NIMBLE models
#'
#' \code{dNmixture_s} and \code{dNmixture_v} provides dynamic occupancy model distributions for NIMBLE models.
#'
#' @name dNmixture
#' @aliases dNmixture_s dNmixture_v rNmixture_s rNmixture_v
#'
#' @author Ben Goldstein, Lauren Ponisio, and Perry de Valpine
#'
#' @param x count observation data, 0 and positive integers
#' @param lambda expected value of the Poisson distribution of true abundance
#' @param prob observation probability for each X
#' @param notZero 0 if datum is structural zero, 1 otherwise. Allows for
#'        zero-inflated Poisson models
#' @param minN the minimum true abundance to sum over. Set to -1 for distribution
#'        to select based on lambda. Ignored if dynamicMinMax = TRUE
#' @param maxN the maximum true abundance to sum over. Set to -1 for distribution
#'        to select based on lambda. Ignored if dynamicMinMax = TRUE
#' @param dynamicMinMax 0 to use input minN and maxN, 1 to ignore
#'        provided minN and maxN and algorithmically select reasonable bounds for N.
#' @param len The length of the x vector. Needed for simulation in \code{rNmixture_*}.
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return probability.
#'        Need not be specified in the model context.
#' @param n number of random draws, each returning a vector of length
#'        \code{len}. Currently only \code{n = 1} is supported, but the
#'        argument exists for standardization of "\code{r}" functions.
#'
#' @details
#' These nimbleFunctions provide distributions that can be used directly in R or
#' in \code{nimble} hierarchical models (via \code{\link[nimble]{nimbleCode}}
#' and \code{\link[nimble]{nimbleModel}}).
#'
#' N-mixture models model abundance at a series of sites over many replicate observations.
#' The likelihood of an observation in site \code{s} at visit \code{t} depends on the
#' true abundance N, which is assumed to be drawn from a Poisson distribution with
#' mean \code{lambda}. It also depends on the probability of detecting an individual
#' (\code{prob} or \code{prob[t]}). The observed count \code{x[t]} has a probability
#' according to the Binomial distribution with size parameter N and probability \code{prob}
#'
#' The distribution has two forms, \code{dNmixture_s} and
#' \code{dNmixture_v}. These differentiate between whether the detection
#' probability \code{prob} is visit-dependent (vector, \code{dNmixture_v})
#' or visit-independent (scalar, dNmixture_s).
#'
#' For more explanation, see
#' \href{../doc/Introduction_to_nimbleEcology.html}{package vignette} (or
#' \code{vignette("Introduction_to_nimbleEcology")}).
#'
#' Compared to writing \code{nimble} models with a discrete latent
#' state of abundance N and a separate scalar datum for each observation time,
#' use of these distributions allows one to directly sum (marginalize) over
#' the discrete latent state N and calculate the probability of all
#' observations for a site jointly.
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
#' minN, maxN, T)}
#'
#' declares that the \code{observedCounts[i, 1:T]} (observed
#' counts for site \code{i}, for example) vector follows an
#' N-mixture model distribution with parameters as indicated,
#' assuming all the parameters have been declared elsewhere in the
#' model. As above, \code{lambda[i]} is the mean of the abundance distribution at site i,
#' \code{prob[i, 1:T]} is a vector of detection probabilities at site i, and
#' \code{T} is the number of observation occasions. This will invoke
#' (something like) the following call to \code{dHMM} when
#' \code{nimble} uses the model such as for MCMC:
#'
#' \code{dNmixture_v(observedCounts[i, 1:T], lambda[i],
#' prob[i, 1:T],
#' minN, maxN, T)}
#'
#' If an algorithm using a \code{nimble} model with this declaration
#' needs to generate a random draw for \code{observedCounts[1:T]}, it
#' will make a similar invocation of \code{rNmixture_v}, with \code{n = 1}.
#'
#' If the observation probabilities are occasion-independent, one would use:
#'
#' \code{observedCounts[1:T] ~ dNmixture_s(observedCounts[i, 1:T], lambda[i],
#' prob[i],
#' minN, maxN, T)}
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
#' @references D. Turek, P. de Valpine and C. J. Paciorek. 2016. Efficient Markov chain Monte
#' Carlo sampling for hierarchical hidden Markov models. Environmental and Ecological Statistics
#' 23:549–564. DOI 10.1007/s10651-016-0353-z
#'
#' @examples
#' \donttest{
#' # Set up constants and initial values for defining the model
#' len <- 5 # length of dataset
#' dat <- c(1,2,0,1,5) # A vector of observations
#' lambda <- 10 # mean abundance
#' prob <- c(0.2, 0.3, 0.2, 0.1, 0.4) # A vector of detection probabilities
#'
#' # Define code for a nimbleModel
#'  nc <- nimbleCode({
#'    x[1:5] ~ dNmixture_v(lambda, prob = [1:5],
#'                         minN = -1, maxN = -1, len = 5)
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
#' nmix_model$calculate()
#' # Use the model for a variety of other purposes...
#' }


NULL
#' @rdname dNmixture
#' @export
dNmixture_v <- nimbleFunction(
    run = function(x = double(1),
                   lambda = double(),
                   prob = double(1),
                   minN = double(),
                   maxN = double(),
                   len = double(),
                   log = integer(0, default = 0)) {
    if (length(x) != len) stop("in dNmixture_v, len must equal length(x).")
    if (length(prob) != len) stop("in dNmixture_v, len must equal length(prob).")

    # Lambda cannot be negative
    if (lambda < 0) {
        if (log) return(-Inf)
        else return(0)
    }

    # Check if there is any data
    # if (is.na.vec(x) | is.nan.vec(x)) {
    #     if (log) return(-Inf)
    #     return(0)
    # }

    ## For each x, the conditional distribution of (N - x | x) is pois(lambda * (1-p))
    ## We determine the lowest N and highest N at extreme quantiles and sum over those.
    if (minN == -1 & maxN == -1) {
      minN <- min(x + qpois(0.00001, lambda * (1 - prob)))
      maxN <- max(x + qpois(0.99999, lambda * (1 - prob)))
      minN <- max( max(x), minN ) ## set minN to at least the largest x
    }

    obsProb <- 0
    if (maxN > minN) { ## should normally be true, but check in case it isn't in some corner case.
    ##    print("counting from ", minN, " to ", maxN, " with lambda = ", lambda)
        for (N in minN:maxN) {
            thisObsProb <- dpois(N, lambda) * prod(dbinom(x, size = N, prob = prob))
            obsProb <- obsProb + thisObsProb
        }
    } else {
        ## return a potentially non-zero obsProb
        N <- max(x)
        obsProb <- dpois(N, lambda) * prod(dbinom(x, size = N, prob = prob))
    }
    if (log) return(log(obsProb))
    else return(obsProb)
    returnType(double())
  })

NULL
#' @rdname dNmixture
#' @export
dNmixture_s <- nimbleFunction(
    run = function(x = double(1),
                   lambda = double(),
                   prob = double(),
                   minN = double(),
                   maxN = double(),
                   len = double(),
                   log = integer(0, default = 0)) {
    if (length(x) != len) stop("in dNmixture_v, len must equal length(x).")

    # Lambda cannot be negative
    if (lambda < 0) {
        if (log) return(-Inf)
        else return(0)
    }

    # Check if there is any data
    # if (is.na.vec(x) | is.nan.vec(x)) {
    #     if (log) return(-Inf)
    #     return(0)
    # }

    ## For each x, the conditional distribution of (N - x | x) is pois(lambda * (1-p))
    ## We determine the lowest N and highest N at extreme quantiles and sum over those.
    if (minN == -1 & maxN == -1) {
      minN <- min(x + qpois(0.00001, lambda * (1 - prob)))
      maxN <- max(x + qpois(0.99999, lambda * (1 - prob)))
      minN <- max( max(x), minN ) ## set minN to at least the largest x
    }

    obsProb <- 0
    if (maxN > minN) { ## should normally be true, but check in case it isn't in some corner case.
    ##    print("counting from ", minN, " to ", maxN, " with lambda = ", lambda)
        for (N in minN:maxN) {
            thisObsProb <- dpois(N, lambda) * prod(dbinom(x, size = N, prob = prob))
            obsProb <- obsProb + thisObsProb
        }
    } else {
        ## return a potentially non-zero obsProb
        N <- max(x)
        obsProb <- dpois(N, lambda) * prod(dbinom(x, size = N, prob = prob))
    }
    if (log) return(log(obsProb))
    else return(obsProb)
    returnType(double())
  })

NULL
#' @rdname dNmixture
#' @export
rNmixture_v <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 prob = double(1),
                 minN = double(),
                 maxN = double(),
                 len = double()) {
    if (n != 1) stop("rNmixture_v only works for n = 1")
    if (length(prob) != len) stop("In rNmixture_v, len must equal length(prob).")
    trueN <- rpois(1, lambda)
    ans <- numeric(len)
    for (i in 1:len) {
      ans[i] <- rbinom(n = 1, size = trueN, prob = prob[i])
    }

    return(ans)
    returnType(double(1))
})

NULL
#' @rdname dNmixture
#' @export
rNmixture_s <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 prob = double(),
                 minN = double(),
                 maxN = double(),
                 len = double()) {
    if (n != 1) stop("rNmixture_v only works for n = 1")
    trueN <- rpois(1, lambda)
    ans <- numeric(len)
    for (i in 1:len) {
      ans[i] <- rbinom(n = 1, size = trueN, prob = prob)
    }

    return(ans)
    returnType(double(1))
})



