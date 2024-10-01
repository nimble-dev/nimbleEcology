
#' Occupancy distribution suitable for use in \code{nimble} models
#'
#' \code{dOcc_*} and \code{rOcc_*} provide occupancy model
#' distributions that can be used directly from R or in \code{nimble}
#' models.
#'
#' @aliases dOcc_s dOcc_v
#'
#' @name dOcc
#'
#' @param x detection/non-detection vector of 0s (not detected) and 1s
#'     (detected).
#' @param probOcc occupancy probability (scalar).
#' @param probDetect detection probability (scalar for \code{dOcc_s},
#'     vector for \code{dOcc_v}).
#' @param len length of detection/non-detection vector (see below).
#' @param log TRUE or 1 to return log probability. FALSE or 0 to
#'     return probability.
#' @param n number of random draws, each returning a vector of length
#'     \code{len}. Currently only \code{n = 1} is supported, but the
#'     argument exists for standardization of "\code{r}" functions.
#'
#' @author Ben Goldstein, Perry de Valpine, and Lauren Ponisio
#'
#' @details
#'
#' These nimbleFunctions provide distributions that can be used directly in R or
#' in \code{nimble} hierarchical models (via \code{\link[nimble]{nimbleCode}}
#' and \code{\link[nimble]{nimbleModel}}).
#'
#' The probability of observation vector \code{x} depends on
#' occupancy probability, \code{probOcc}, and detection probability,
#' \code{probDetect} or \code{probDetect[1:T]}.
#'
#' The letter following the 'dOcc_' indicates whether detection probability is
#' scalar (s, meaning \code{probDetect} is detection probability for every
#' \code{x[t]}) or vector (v, meaning \code{probDetect[t]} is detection
#' probability for \code{x[t]}).
#'
#' When used directly from R, the \code{len} argument to \code{dOcc_*} is not
#' necessary. It will default to the length of \code{x}.  When used in
#' \code{nimble} model code (via \code{nimbleCode}), \code{len} must be provided
#' (even though it may seem redundant).
#'
#' For more explanation, see package vignette
#' (\code{vignette("Introduction_to_nimbleEcology")}).
#'
#' Compared to writing \code{nimble} models with a discrete latent state for
#' true occupancy status and a separate scalar datum for each observation, use
#' of these distributions allows one to directly sum (marginalize) over the
#' discrete latent state and calculate the probability of all observations from
#' one site jointly.
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
#' \code{detections[i, 1:T] ~ dOcc_s(occupancyProbability,
#' detectionProbability, T)}
#'
#' declares that \code{detections[i, 1:T]} (detection history at site \code{i},
#' for example) follows an occupancy distribution with parameters as indicated,
#' assuming all the parameters have been declared elsewhere in the model.  This
#' will invoke (something like) the following call to \code{dOcc_s} when
#' \code{nimble} uses the model such as for MCMC:
#'
#' \code{dOcc_s(detections[i, 1:T], occupancyProbability,
#' detectionProbability, len = T, log = TRUE)}
#'
#' If an algorithm using a \code{nimble} model with this declaration
#' needs to generate a random draw for \code{detections[i, 1:T]}, it
#' will make a similar invocation of \code{rOcc_s}, with \code{n = 1}.
#'
#' If the detection probabilities are time-dependent, use:
#'
#' \code{detections[i, 1:T] ~ dOcc_v(occupancyProbability,
#' detectionProbability[1:T], len = T)}
#'
#' @section Notes for use with automatic differentiation:
#'
#' The \code{dOcc_*} distributions should all work for models and algorithms
#' that use nimble's automatic differentiation (AD) system. In that system, some
#' kinds of values are "baked in" (cannot be changed) to the AD calculations
#' from the first call, unless and until the AD calculations are reset. For the
#' \code{dOcc_*} distributions, the lengths of vector inputs are baked in. These
#' can be different for different iterations through a for loop (or nimble model
#' declarations with different indices, for example), but the lengths for each
#' specific iteration will be "baked in" after the first call. \bold{It is
#' safest if one can assume that \code{x} are data and are not going to change.}
#'
#' @return
#'
#' For \code{dOcc_*}: the probability (or likelihood) or log probability of observation vector \code{x}.
#'
#' For \code{rOcc_*}: a simulated detection history, \code{x}.
#'
#' @seealso For dynamic occupancy models, see documentation for
#'   \code{\link{dDynOcc}}.
#' @examples
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
    logProb_x_given_occupied <- 0
    prob_x_given_unoccupied <- 1
    for(i in 1:length(x)) {
      xi <- ADbreak(x[i])
      if(!is.na(xi)) { # Handle missing values
        logProb_x_given_occupied <- logProb_x_given_occupied + dbinom(x[i], prob = probDetect, size = 1, log = TRUE)
        if(xi==1) prob_x_given_unoccupied <- 0
      }
    }
    prob_x <- exp(logProb_x_given_occupied) * probOcc + prob_x_given_unoccupied * (1 - probOcc)
    if (log) return(log(prob_x))
    return(prob_x)
  }, buildDerivs = list(run = list(ignore = c("i", "xi")))
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
    logProb_x_given_occupied <- 0
    prob_x_given_unoccupied <- 1
    for(i in 1:length(x)) {
      xi <- ADbreak(x[i])
      if(!is.na(xi)) { # Handle missing values
        logProb_x_given_occupied <- logProb_x_given_occupied + dbinom(x[i], prob = probDetect[i], size = 1, log = TRUE)
        if(xi==1) prob_x_given_unoccupied <- 0
      }
    }
    prob_x <- exp(logProb_x_given_occupied) * probOcc + prob_x_given_unoccupied * (1 - probOcc)
    if (log) return(log(prob_x))
    return(prob_x)
  }, buildDerivs = list(run = list(ignore = c("i", "xi")))
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
