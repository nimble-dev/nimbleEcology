#' Occupancy distribution for use in NIMBLE models
#'
#' \code{dOcc_s} and \code{dOcc_v} provide occupancy model distributions for NIMBLE models.
#' The 's' version uses time-independent (scalar) detection probability while the 'v' version
#' uses time-dependent (vector) detection probability.
#' Compared to writing NIMBLE models with a discrete latent state for true occupancy status and
#' a separate scalar datum for each observation,
#' use of these distributions allows
#' one to directly sum over the discrete latent state and calculate the probability of
#' all observations from one site jointly.
#'
#' @aliases dOcc rOcc rOcc_s dOcc_v rOcc_v
#'
#' @export
#'
#' @param x detection/non-detection vector of 0s (not detected) and 1s (detected).
#' @param probOcc occupancy probability (scalar).
#' @param probDetect detection probability (scalar for \code{dOcc_s}, vector for \code{dOcc_v}).
#' @param l length of detection/non-detection vector (ignored for "d" functions, needed for "r" functions).
#' @param log TRUE (return log probability) or FALSE (return probability)
#'
#' @author Perry de Valpine
#'
#' @details These nimbleFunctions provide distributions that can be used in code (via \link{nimbleCode})
#' for \link{nimbleModel}.
#'
#' These are written in the format of user-defined distributions to extend NIMBLE's
#' use of the BUGS model language.  More information about writing user-defined distributions can be found
#' in the NIMBLE User Manual at \code{https://r-nimble.org}.
#'
#' The first argument to a "d" function is always named \code{x} and is given on the
#' left-hand side of a (stochastic) model declaration in the BUGS model language (used by NIMBLE).
#' When using these distributions in a NIMBLE model, the user
#' should not provide the \code{log} argument.  (It is always set to \code{TRUE} when used
#' in a NIMBLE model.)
#'
#' For example, in a NIMBLE model,
#'
#' \code{detections[1:T] ~ dOcc_s(occupancyProbability, detectionProbability)}
#'
#' declares that the \code{detections[1:T]} vector follows an occupancy model distribution
#' with parameters as indicated, assuming all the parameters have been declared elsewhere in the model.
#'
#' If the detection probabilities are time-dependent, one would use:
#'
#'\code{detections[1:T] ~ dOcc_s(occupancyProbability, detectionProbability[1:T])}
#'
#' @seealso For dynamic occupancy models, see \link{dDynOcc}.
dOcc_s <- nimbleFunction(
  run = function(x = double(1),
                 probOcc = double(),
                 probDetect = double(),
                 log = logical(0, default = 0)) {
    returnType(double(0))
    logProb_x_given_occupied <- sum(dbinom(x, prob = probDetect, size = 1, log = TRUE))
    prob_x_given_unoccupied <- sum(x) == 0
    prob_x <- exp(logProb_x_given_occupied) * probOcc + prob_x_given_unoccupied * (1 - probOcc)
    if (log) return(log(prob_x))
    return(prob_x)
  }
)

dOcc_v <- nimbleFunction(
  run = function(x = double(1),
                 probOcc = double(0),
                 probDetect = double(1),
                 log = logical(0, default = 0)) {
    returnType(double(0))
    logProb_x_given_occupied <- sum(dbinom(x, prob = probDetect, size = 1, log = TRUE))
    prob_x_given_unoccupied <- sum(x) == 0
    prob_x <- exp(logProb_x_given_occupied) * probOcc + prob_x_given_unoccupied * (1 - probOcc)
    if (log) return(log(prob_x))
    return(prob_x)
  }
)

rOcc_s <- nimbleFunction(
  run = function(n = integer(),
                 probOcc = double(),
                 probDetect = double()) {
    returnType(double(1))
    k <- length(probDetect)
    z <- rbinom(1, prob = probOcc, size = 1)
    if(z == 0) return(numeric(k))
    return(rbinom(k, prob = probDetect, size = 1))
  }
)

rOcc_v <- nimbleFunction(
  run = function(n = integer(),
                 probOcc = double(0),
                 probDetect = double(1)) {
    returnType(double(1))
    k <- length(probDetect)
    z <- rbinom(1, prob = probOcc, size = 1)
    if (z == 0) return(numeric(k))
    return(rbinom(k, prob = probDetect, size = 1))
  }
)


# registerDistributions(list(
#   dOcc_s = list(
#     BUGSdist = "dOcc_s(probOcc, probDetect)",
#     Rdist = "dOcc_s(probOcc, probDetect)",
#     discrete = TRUE,
#     types = c('value = double(1)', 'probOcc = double(0)', 'probDetect = double(0)'),
#     pqAvail = FALSE))
# )
