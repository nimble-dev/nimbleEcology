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
#' @param probOcc occupancy probability
#' @param probDetect detection probability
#' @param l length of detection/non-detection history (needed for rOcc)
#' @param log TRUE (return log probability) or FALSE (return probability)
#'
#' @details To be filled in.
dOcc_s <- nimbleFunction(
  run = function(x = double(1),
                 probOcc = double(0),
                 probDetect = double(0),
                 log = logical(0, default = 0)) {
    returnType(double(0))
    logProb_x_given_occupied <- sum(dbinom(x, prob = probDetect, size = 1, log = TRUE))
    prob_x_given_unoccupied <- sum(x) == 0
    prob_x <- exp(logProb_x_given_occupied) * probOcc + prob_x_given_unoccupied * (1-probOcc)
    if(log) return(log(prob_x))
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
    prob_x <- exp(logProb_x_given_occupied) * probOcc + prob_x_given_unoccupied * (1-probOcc)
    if(log) return(log(prob_x))
    return(prob_x)
  }
)

rOcc_s <- nimbleFunction(
  run = function(n = integer(),
                 probOcc = double(0),
                 probDetect = double(0)) {
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
    if(z == 0) return(numeric(k))
    return(rbinom(k, prob = probDetect, size = 1))
  }
)
