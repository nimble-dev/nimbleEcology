# dNmixture
#' N-mixture distribution for use in NIMBLE models
#'
#' \code{dNmixture} provides dynamic occupancy model distributions for NIMBLE models.
#'
#' N-mixture models
#' model abundance at a series of sites over many replicate observations. The likelihood of an observation in site s
#' at time t depends on the
#'
#' Compared to writing NIMBLE models with a discrete latent state for true occupancy status and
#' a separate scalar datum for each observation,
#' use of these distributions allows
#' one to directly sum over the discrete latent state and calculate the probability of
#' all observations from one site jointly.
#'
#' @name dNmixture
#' @aliases dDynOcc_ss dDynOcc_sv dDynOcc_vs dDynOcc_vv
#'
#' @param x count observation data, 0 and positive integers
#' @param lambda expected value of the Poisson distribution of true abundance
#' @param prob observation probability for each X
#' @param notZero 0 if datum is structural zero, 1 otherwise. Allows for zero-inflated Poisson models
#' @param minN the minimum true abundance to sum over. Set to -1 for distribution to select based on lambda
#' @param maxN the maximum true abundance to sum over. Set to -1 for distribution to select based on lambda
#'
#' @author Lauren Ponisio, Ben Goldstein, Perry de Valpine
#'
#' These are written in the format of user-defined distributions to extend NIMBLE's
#' use of the BUGS model language. More information about writing user-defined distributions can be found
#' in the NIMBLE User Manual at \code{https://r-nimble.org}.
#'
#' The first argument to a "d" function is always named \code{x} and is given on the
#' left-hand side of a (stochastic) model declaration in the BUGS model language (used by NIMBLE).
#' When using these distributions in a NIMBLE model, the user
#' should not provide the \code{log} argument. (It is always set to \code{TRUE} when used
#' in a NIMBLE model.)
#'
#' For example, in a NIMBLE model,
#'
#' \code{detections[1:T] ~ dNmixture(lambda = lambda, prob = p[1:T], notZero = 1, minN = -1, maxN = -1)}
#'
#' declares that the \code{detections[1:T]} vector follows an N-mixture model distribution
#' with parameters as indicated, assuming all the parameters have been declared elsewhere in the model.
#'
#'
#' @seealso For regular occupancy models, see documentation for dOcc.
#' @examples
#' \dontrun{}


NULL
#' @rdname dNmixture
#' @export
dNmixture <- nimbleFunction(
    run = function(x = double(1),
                   lambda = double(),
                   prob = double(1),
                   notZero = double(),
                   minN = double(),
                   maxN = double(),
                   log = integer(0, default = 0)) {
    if (lambda < 0) {
        if (log) return(-Inf)
        else return(0)
    }
    if (is.na.vec(x) | is.nan.vec(x)) {
        if (log) return(-Inf)
        return(0)
    }
    if (notZero == 0) { ## It is structural zero
        if (all(x == 0)) {
            if (log) return(0)
            else return(1)
        } else {
            if (log) return(-Inf)
            else return(0)
        }
    }
    ## It is not a structural zero.
    ##
    ## For each x, the conditional distribution of (N - x | x) is pois(lambda * (1-p))
    ## We determine the lowest N and highest N at extreme quantiles and sum over those.
    minN <- min(x + qpois(0.00001, lambda * (1 - prob)))
    maxN <- max(x + qpois(0.99999, lambda * (1 - prob)))
    minN <- max( max(x), minN ) ## set minN to at least the largest x

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
    returnType(double(0))
  })
