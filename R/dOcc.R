#' Occupancy distribution for use in NIMBLE models
#'
#' \code{dOcc_**} provides occupancy model distributions for NIMBLE models.
#'
#' Likelihood of observation x[t] depends on an occupancy probability \code{probOcc[t]} and
#' detection probability \code{probDetect[t]}.
#' The pair of letters following the 'dOcc_' indicates whether the occupancy probability
#' and detection probability are scalar (s) or vector (v). For example, dOcc_sc takes scalar
#' occupancy probability with a vector of detection probabilities.
#'
#' Compared to writing NIMBLE models with a discrete latent state for true occupancy status and
#' a separate scalar datum for each observation,
#' use of these distributions allows
#' one to directly sum over the discrete latent state and calculate the probability of
#' all observations from one site jointly.
#'
#' @aliases dOcc_s dOcc_v
#'
#' @name dOcc
#'
#' @param x detection/non-detection vector of 0s (not detected) and 1s (detected).
#' @param probOcc occupancy probability (scalar for \code{dOcc_s*}, vector for \code{dOcc_v*}).
#' @param probDetect detection probability (scalar for \code{dOcc_*s}, vector for \code{dOcc_*v}).
#' @param len length of detection/non-detection vector (ignored for "d" functions, needed for "r" functions).
#' @param log TRUE (return log probability) or FALSE (return probability)
#' @param n length of random sequence
#'
#' @author Ben Goldstein and Perry de Valpine
#'
#' @details These nimbleFunctions provide distributions that can be used in code (via \link{nimbleCode})
#' for \link{nimbleModel}.
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
#' \code{detections[1:T] ~ dOcc_s(occupancyProbability, detectionProbability, len = T)}
#'
#' declares that the \code{detections[1:T]} vector follows an occupancy model distribution
#' with parameters as indicated, assuming all the parameters have been declared elsewhere in the model.
#'
#' If the detection probabilities are time-dependent, one would use:
#'
#' \code{detections[1:T] ~ dOcc_sv(occupancyProbability, detectionProbability[1:T], len = T)}
#'
#' @seealso For dynamic occupancy models, see documentation for \link{dDynOcc}.
NULL
#' @examples
#' \dontrun{
#' # Set up constants and initial values for defining the model
#' dat <- c(1,1,0,0) # A vector of observations
#' probOcc <- 0.6
#' probDetect <- 0.4
#'
#'
#' # Define code for a nimbleModel
#' nc <- nimbleCode({
#'   x[1:4] ~ dOcc_s(probOcc, probDetect, len = 4)
#'   probOcc ~ dunif(0,1)
#'   probDetect ~ dunif(0,1)
#' })
#'
#' # Build the model, providing data and initial values
#' Occ_model <- nimbleModel(nc, data = list(x = dat),
#'                          inits = list(probOcc = probOcc,
#'                                       probDetect = probDetect))
#'
#' # Calculate log probability of data from the model
#' Occ_model$calculate()
#' # Use the model for a variety of other purposes...
#' }
#' @export
#' @rdname dOcc
dOcc_s <- nimbleFunction(
  run = function(x = double(1),
                 probOcc = double(0),
                 probDetect = double(0),
                 len = integer(0, default = 0),
                 log = logical(0, default = 0)) {
    if (len != 0) if (len != length(x)) stop("Argument 'len' must match length of data, or be 0.")
    returnType(double(0))
    logProb_x_given_occupied <- sum(dbinom(x, prob = probDetect, size = 1, log = TRUE))
    prob_x_given_unoccupied <- sum(x) == 0
    prob_x <- exp(logProb_x_given_occupied) * probOcc + prob_x_given_unoccupied * (1 - probOcc)
    if (log) return(log(prob_x))
    return(prob_x)
  }
)

#' @export
#' @rdname dOcc
dOcc_v <- nimbleFunction(
  run = function(x = double(1),
                 probOcc = double(0),
                 probDetect = double(1),
                 len = integer(0, default = 0),
                 log = logical(0, default = 0)) {
    if (len != 0) if (len != length(x)) stop("Argument 'len' must match length of data, or be 0.")
    if (length(x) != length(probDetect)) stop("Length of data does not match length of detection vector.")
    returnType(double(0))
    logProb_x_given_occupied <- sum(dbinom(x, prob = probDetect, size = 1, log = TRUE))
    prob_x_given_unoccupied <- sum(x) == 0
    prob_x <- exp(logProb_x_given_occupied) * probOcc + prob_x_given_unoccupied * (1 - probOcc)
    if (log) return(log(prob_x))
    return(prob_x)
  }
)

#' @export
#' @rdname dOcc

rOcc_s <- nimbleFunction(
  run = function(n = integer(),
                 probOcc = double(0),
                 probDetect = double(0),
                 len = integer(0, default = 0)) {
    if (len == 0) stop("Argument 'len' must be given for rOcc_s (or if a nimble model with dOcc_s is used for simulation).")
    returnType(double(1))
    u <- runif(1, 0, 1)
    if (u > probOcc) return(numeric(0, length = len))
    return(rbinom(len, prob = probDetect, size = 1))
  }
)

#' @export
#' @rdname dOcc
rOcc_v <- nimbleFunction(
  run = function(n = integer(),
                 probOcc = double(0),
                 probDetect = double(1),
                 len = integer(0, default = 0)) {
    if(len != 0) {
      if (len != length(probDetect)) {
        stop("If argument 'len' is given, it must match length of probDetect.")
      }
    }
    returnType(double(1))
    k <- length(probDetect)
    u <- runif(1, 0, 1)
    if(u > probOcc) return(numeric(0, length = k))
    return(rbinom(k, prob = probDetect, size = 1))
  }
)

registerDistributions(list(
  dOcc_s = list(
    BUGSdist = "dOcc_s(probOcc, probDetect, len)",
    Rdist = "dOcc_s(probOcc, probDetect, len)",
    discrete = TRUE,
    types = c('value = double(1)', 'probOcc = double(0)', 'probDetect = double(0)', 'len = integer(0)'),
    pqAvail = FALSE)))

registerDistributions(list(
  dOcc_v = list(
    BUGSdist = "dOcc_v(probOcc, probDetect, len)",
    Rdist = c("dOcc_v(probOcc, probDetect, len)"),
    discrete = TRUE,
    types = c('value = double(1)', 'probOcc = double(0)', 'probDetect = double(1)', 'len = integer(0)'),
    pqAvail = FALSE)))
