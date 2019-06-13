#' Occupancy distribution for use in NIMBLE models
#'
#' \code{dOcc_**} provides occupancy model distributions for NIMBLE models.
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
#' @aliases dOcc_ss dOcc_sv dOcc_vs dOcc_vv rOcc_ss rOcc_sv rOcc_vs rOcc_vv
#'
#' @name dOcc
#'
#' @export
#'
#' @param x detection/non-detection vector of 0s (not detected) and 1s (detected).
#' @param probOcc occupancy probability (scalar for \code{dOcc_\*s}, vector for \code{dOcc_\*v}).
#' @param probDetect detection probability (scalar for \code{dOcc_\*s}, vector for \code{dOcc_\*v}).
#' @param len length of detection/non-detection vector (ignored for "d" functions, needed for "r" functions).
#' @param log TRUE (return log probability) or FALSE (return probability)
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
#' \code{detections[1:T] ~ dOcc_ss(occupancyProbability, detectionProbability)}
#'
#' declares that the \code{detections[1:T]} vector follows an occupancy model distribution
#' with parameters as indicated, assuming all the parameters have been declared elsewhere in the model.
#'
#' If the detection probabilities are time-dependent, one would use:
#'
#' \code{detections[1:T] ~ dOcc_sv(occupancyProbability, detectionProbability[1:T])}
#'
#' @seealso For dynamic occupancy models, see \link{dDynOcc}.

#' @export
#' @rdname dOcc
dOcc_ss <- nimbleFunction(
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
dOcc_sv <- nimbleFunction(
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
dOcc_vs <- nimbleFunction(
  run = function(x = double(1),
                 probOcc = double(1),
                 probDetect = double(0),
                 len = integer(0, default = 0),
                 log = logical(0, default = 0)) {
    if (len != 0) if (len != length(x)) stop("Argument 'len' must match length of data, or be 0.")
    if (length(x) != length(probOcc)) stop("Length of data does not match length of occupancy vector.")

    returnType(double(0))
    prob_x_given_occupied <- dbinom(x, prob = probDetect, size = 1, log = FALSE)
    prob_x_given_unoccupied <- x == 0
    log_prob_x <- sum(log(prob_x_given_occupied * probOcc + prob_x_given_unoccupied * (1 - probOcc)))

    if (log) return(log_prob_x)
    return(exp(log_prob_x))
  }
)

#' @export
#' @rdname dOcc
dOcc_vv <- nimbleFunction(
  run = function(x = double(1),
                 probOcc = double(1),
                 probDetect = double(1),
                 len = integer(0, default = 0),
                 log = logical(0, default = 0)) {
    if (len != 0) if (len != length(x)) stop("Argument 'len' must match length of data, or be 0.")
    if (length(x) != length(probOcc)) stop("Length of data does not match length of occupancy vector.")
    if (length(x) != length(probDetect)) stop("Length of data does not match length of detection vector.")

    returnType(double(0))
    prob_x_given_occupied <- dbinom(x, prob = probDetect, size = 1, log = FALSE)
    prob_x_given_unoccupied <- x == 0
    log_prob_x <- sum(log(prob_x_given_occupied * probOcc + prob_x_given_unoccupied * (1 - probOcc)))
    if (log) return(log_prob_x)
    return(exp(log_prob_x))
  }
)

#' @export
#' @rdname dOcc
rOcc_ss <- nimbleFunction(
  run = function(n = integer(),
                 probOcc = double(0),
                 probDetect = double(0),
                 len = integer(0, default = 0)) {
    if (len == 0) stop("Argument 'len' must be specified in rOcc_ss form.")
    returnType(double(1))
    k <- len
    z <- rbinom(1, prob = probOcc, size = 1)
    if (z == 0) return(numeric(k))
    return(rbinom(k, prob = probDetect, size = 1))
  }
)

#' @export
#' @rdname dOcc
rOcc_sv <- nimbleFunction(
  run = function(n = integer(),
                 probOcc = double(0),
                 probDetect = double(1),
                 len = integer(0, default = 0)) {
    if (len != length(probDetect)) {
      stop("If specified, argument 'len' must match length of probDetect.")
    }
    returnType(double(1))
    k <- length(probDetect)
    z <- rbinom(1, prob = probOcc, size = 1)
    if (z == 0) return(numeric(k))
    return(rbinom(k, prob = probDetect, size = 1))
  }
)

#' @export
#' @rdname dOcc
rOcc_vs <- nimbleFunction(
  run = function(n = integer(0),
                 probOcc = double(1),
                 probDetect = double(0),
                 len = integer(0, default = 0)) {
    if (len != length(probOcc)) {
      stop("If specified, argument 'len' must match length of probOcc")
    }
    returnType(double(1))
    k <- length(probOcc)
    z <- rbinom(k, prob = probOcc, size = 1)
    return(rbinom(k, prob = probDetect, size = 1) * z)
  }
)

#' @export
#' @rdname dOcc
rOcc_vv <- nimbleFunction(
  run = function(n = integer(),
                 probOcc = double(1),
                 probDetect = double(1),
                 len = integer(0, default = 0)) {
    if (length(probDetect) != length(probOcc)) {
      stop("Occupancy and detection vectors are of different lengths.")
    }
    if (len != length(probDetect)) {
      stop("If specified, argument 'len' must match length of probDetect and probOcc.")
    }
    returnType(double(1))
    k <- length(probDetect)
    z <- rbinom(k, prob = probOcc, size = 1)
    return(rbinom(k, prob = probDetect, size = 1) * z)
  }
)


registerDistributions(list(
  dOcc_ss = list(
    BUGSdist = "dOcc_ss(probOcc, probDetect, len)",
    Rdist = "dOcc_ss(probOcc, probDetect, len)",
    discrete = TRUE,
    types = c('value = double(1)', 'probOcc = double(0)', 'probDetect = double(0)', 'len = integer(0)'),
    pqAvail = FALSE)))

registerDistributions(list(
  dOcc_sv = list(
    BUGSdist = "dOcc_sv(probOcc, probDetect, len)",
    Rdist = "dOcc_sv(probOcc, probDetect, len)",
    discrete = TRUE,
    types = c('value = double(1)', 'probOcc = double(0)', 'probDetect = double(1)', 'len = integer(0)'),
    pqAvail = FALSE)))

registerDistributions(list(
  dOcc_vs = list(
    BUGSdist = "dOcc_vs(probOcc, probDetect, len)",
    Rdist = "dOcc_vs(probOcc, probDetect, len)",
    discrete = TRUE,
    types = c('value = double(1)', 'probOcc = double(1)', 'probDetect = double(0)', 'len = integer(0)'),
    pqAvail = FALSE)))

registerDistributions(list(
  dOcc_vv = list(
    BUGSdist = "dOcc_vv(probOcc, probDetect, len)",
    Rdist = "dOcc_vv(probOcc, probDetect, len)",
    discrete = TRUE,
    types = c('value = double(1)', 'probOcc = double(1)', 'probDetect = double(1)', 'len = integer(0)'),
    pqAvail = FALSE)))
