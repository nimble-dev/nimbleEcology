#' Dynamic Hidden Markov Model distribution for use in \code{nimble} models
#'
#' \code{dDHMM} and \code{dDHMMo} provide Dynamic hidden Markov model
#' distributions for \code{nimble} models.
#'
#' @name dDHMM
#' @aliases dDHMM dDHMMo rDHMM rDHMMo
#' @author Perry de Valpine, Daniel Turek, and Ben Goldstein
#' @export
#'
#' @param x vector of observations, each one a positive integer
#'     corresponding to an observation state
#'     (one value of which could can correspond to "not observed", and
#'     another value of which can correspond to "dead" or
#'     "removed from system").
#' @param init vector of initial state probabilities. Must sum to 1
#' @param probObs time-independent matrix (\code{dDHMM} and
#'     \code{rHMM}) or time-dependent 3D array (\code{dDHMMo} and
#'     \code{rHMMo}) of observation probabilities.
#'     First two dimensions of \code{probObs} are of size x (number of possible
#'     system states) x (number of possible observation classes). \code{dDHMMo}
#'     and \code{rHMMo} expects an additional third dimension of size (number of
#'     observation times)
#' @param probTrans time-dependent array of system state transition
#' probabilities. Dimension of \code{probTrans} is (number of possible
#' system states) x  (number of possible system states)
#' x (number of observation times)
#' @param len length of observations (needed for rDHMM)
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return probability
#' @param n number of random draws, each returning a vector of length
#'     \code{len}. Currently only \code{n = 1} is supported, but the
#'     argument exists for standardization of "\code{r}" functions
#'
#' @details
#' These nimbleFunctions provide distributions that can be used directly in R or
#' in \code{nimble} hierarchical models (via \code{\link[nimble]{nimbleCode}}
#' and \code{\link[nimble]{nimbleModel}}).
#'
#' The probability (or likelihood) of observation \code{x[t, o]} depends on
#' the previous true latent state, the time-dependent probability of
#' transitioning to a new state \code{probTrans}, and the probability of
#' observation states given the true latent state \code{probObs}.
#'
#' The distribution has two forms, \code{dDHMM} and \code{dDHMMo}. \code{dDHMM}
#' takes a time-independent observation probability matrix with dimension
#' S x O, while \code{dDHMMo} expects a three-dimensional array of time-dependent
#' observation probabilities with dimension S x O x T, where O is the number of
#' possible occupancy states, S is the number of true latent states, and T is
#' the number of time intervals.
#'
#' \code{probTrans} has dimension S x S x (T - 1). \code{probTrans}[i, j, t] is
#' the probability that an individual in state \code{i} at time \code{t} takes on
#' state \code{j} at time \code{t+1}. The length of the third dimension may be greater
#' than (T - 1) but all values indexed greater than T - 1 will be ignored.
#'
#' \code{initStates} has length S. \code{initStates[i]} is the
#' probability of being in state \code{i} at the first observation time.
#'
#' For more explanation, see
#' \href{../doc/Introduction_to_nimbleEcology.html}{package vignette} (or
#' \code{vignette("Introduction_to_nimbleEcology")}).
#'
#' Compared to writing \code{nimble} models with a discrete true latent state
#' and a separate scalar datum for each observation, use
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
#' For example, in a NIMBLE model,
#'
#' \code{observedStates[1:T] ~ dDHMM(initStates[1:S],
#' observationProbs[1:S, 1:O],
#' transitionProbs[1:S, 1:S, 1:(T-1)], T)}
#'
#' declares that the \code{observedStates[1:T]} vector follows a dynamic hidden
#' Markov model distribution with parameters as indicated, assuming all the
#' parameters have been declared elsewhere in the model. In this case, \code{S}
#' is the number of system states, \code{O} is the number of observation
#' classes, and \code{T} is the number of observation occasions.This
#' will invoke (something like) the following call to \code{dDHMM} when
#' \code{nimble} uses the model such as for MCMC:
#'
#' \code{rDHMM(observedStates[1:T], initStates[1:S],
#' observationProbs[1:S, 1:O],
#' transitionProbs[1:S, 1:S, 1:(T-1)], T, log = TRUE)}
#'
#' If an algorithm using a \code{nimble} model with this declaration
#' needs to generate a random draw for \code{observedStates[1:T]}, it
#' will make a similar invocation of \code{rDHMM}, with \code{n = 1}.
#'
#' If the observation probabilities are time-dependent, one would use:
#'
#' \code{observedStates[1:T] ~
#' dDHMMo(initStates[1:S], observationProbs[1:S, 1:O, 1:(T-1)],
#' transitionProbs[1:S, 1:S, 1:(T-1)], T)}
#'
#' @return
#' For \code{dDHMM} and \code{dDHMMo}: the probability (or likelihood) or log
#' probability of observation vector \code{x}.
#' For \code{rDHMM} and \code{rDHMMo}: a simulated detection history, \code{x}.

#' @seealso For hidden Markov models with time-independent transitions,
#' see \link{dHMM} and \link{dHMMo}.
#' For simple capture-recapture, see \link{dCJS}.
#'
#' @references D. Turek, P. de Valpine and C. J. Paciorek. 2016. Efficient Markov chain Monte
#' Carlo sampling for hierarchical hidden Markov models. Environmental and Ecological Statistics
#' 23:549–564. DOI 10.1007/s10651-016-0353-z
#' @examples
#' \donttest{
#' # Set up constants and initial values for defining the model
#' dat <- c(1,2,1,1) # A vector of observations
#' init <- c(0.4, 0.2, 0.4) # A vector of initial state probabilities
#' probObs <- t(array( # A matrix of observation probabilities
#'        c(1, 0,
#'          0, 1,
#'          0.8, 0.2), c(2, 3)))
#'
#' probTrans <- array(rep(0.5, 27), # A matrix of time-indexed transition probabilities
#'             c(3,3,3))
#'
#' # Define code for a nimbleModel
#'  nc <- nimbleCode({
#'    x[1:4] ~ dDHMM(init[1:3], probObs = probObs[1:3, 1:2],
#'                   probTrans = probTrans[1:3, 1:3, 1:3], len = 4)
#'
#'    for (i in 1:3) {
#'      init[i] ~ dunif(0,1)
#'
#'      for (j in 1:3) {
#'        for (t in 1:3) {
#'          probTrans[i,j,t] ~ dunif(0,1)
#'        }
#'      }
#'
#'      probObs[i, 1] ~ dunif(0,1)
#'      probObs[i, 2] <- 1 - probObs[i,1]
#'    }
#'  })
#'
#' # Build the model, providing data and initial values
#' DHMM_model <- nimbleModel(nc,
#'                           data = list(x = dat),
#'                           inits = list(init = init,
#'                                        probObs = probObs,
#'                                        probTrans = probTrans))
#' # Calculate log probability of x from the model
#' DHMM_model$calculate()
#' # Use the model for a variety of other purposes...
#' }

NULL

#' @export
#' @rdname dDHMM
dDHMM <- nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 init = double(1),
                 probObs = double(2),
                 probTrans = double(3),
                 len = double(),## length of x (needed as a separate param for rDHMM)
                 log = integer(0, default = 0)) {
    if (length(init) != dim(probObs)[1]) stop("Length of init does not match nrow of probObs in dDHMM.")
    if (length(init) != dim(probTrans)[1]) stop("Length of init does not match dim(probTrans)[1] in dDHMM.")
    if (length(init) != dim(probTrans)[2]) stop("Length of init does not match dim(probTrans)[2] in dDHMM.")
    if (length(x) != len) stop("Length of x does not match len in dDHMM.")
    if (len - 1 != dim(probTrans)[3]) stop("len - 1 does not match dim(probTrans)[3] in dDHMM.")

    pi <- init # State probabilities at time t=1
    logL <- 0
    nObsClasses <- dim(probObs)[2]
    lengthX <- length(x)
    for (t in 1:lengthX) {
      if (x[t] > nObsClasses) stop("Invalid value of x[t] in dDHMM.")
      Zpi <- probObs[, x[t]] * pi # Vector of P(state) * P(observation class x[t] | state)
      sumZpi <- sum(Zpi)    # Total P(observed as class x[t])
      logL <- logL + log(sumZpi)  # Accumulate log probabilities through time
      if (t != lengthX) pi <- (probTrans[,,t] %*% asCol(Zpi) / sumZpi)[ ,1] # State probabilities at t+1
    }

    returnType(double())
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @export
#' @rdname dDHMM
dDHMMo <- nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 init = double(1),##
                 probObs = double(3),
                 probTrans = double(3),
                 len = double(),## length of x (needed as a separate param for rDHMM)
                 log = integer(0, default = 0)) {
    if (length(init) != dim(probObs)[1]) stop("Length of init does not match ncol of probObs in dDHMMo.")
    if (length(init) != dim(probTrans)[1]) stop("Length of init does not match dim(probTrans)[1] in dDHMMo.")
    if (length(init) != dim(probTrans)[2]) stop("Length of init does not match dim(probTrans)[2] in dDHMMo.")
    if (length(x) != len) stop("Length of x does not match len in dDHMM.")
    if (len - 1 > dim(probTrans)[3]) stop("dim(probTrans)[3] does not match len - 1 in dDHMMo.")
    if (len != dim(probObs)[3]) stop("dim(probObs)[3] does not match len in dDHMMo.")

    pi <- init # State probabilities at time t=1
    logL <- 0
    nObsClasses <- dim(probObs)[2]
    lengthX <- length(x)
    for (t in 1:lengthX) {
      if (x[t] > nObsClasses) stop("Invalid value of x[t] in dDHMM.")
      Zpi <- probObs[, x[t], t] * pi # Vector of P(state) * P(observation class x[t] | state)
      sumZpi <- sum(Zpi)    # Total P(observed as class x[t])
      logL <- logL + log(sumZpi)  # Accumulate log probabilities through time
      if (t != lengthX)   pi <- (probTrans[,,t] %*% asCol(Zpi) / sumZpi)[ ,1] # State probabilities at t+1
    }
    returnType(double())
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @export
#' @rdname dDHMM
rDHMM <- nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 init = double(1),
                 probObs = double(2),
                 probTrans = double(3),
                 len = double()) {
    nStates <- length(init)
    if (nStates != dim(probObs)[1]) stop("Length of init does not match nrow of probObs in dDHMM.")
    if (nStates != dim(probTrans)[1]) stop("Length of init does not match dim(probTrans)[1] in dDHMM.")
    if (nStates != dim(probTrans)[2]) stop("Length of init does not match dim(probTrans)[2] in dDHMM.")
    if (len - 1 > dim(probTrans)[3]) stop("len - 1 does not match dim(probTrans)[3] in dDHMM.")

    returnType(double(1))
    ans <- numeric(len)

    trueState <- rcat(1, init)
    for (i in 1:len) {
      # Detect based on the true state
      ans[i] <- rcat(1, probObs[trueState,])
      # Transition to a new true state
      if (i != len) {
        trueState <- rcat(1, probTrans[trueState, , i])
    }
  }
  return(ans)
})
# rDHMM <- nimbleFunction(
#   run = function(n = integer(),    ## Observed capture (state) history
#                  init = double(1), ## probabilities of state at time 1
#                  probObs = double(2),
#                  probTrans = double(3),
#                  len = integer()) {
#   returnType(double(1))
#   ans <- numeric(len)
#
#   trueState <- rmulti(n = 1, size = 1, prob = init)
#
#   for (t in 1:len) {
#     thisstate <- which(trueState == 1)
#     ans[t] <- which(rmulti(1, size = 1, prob = probObs[, thisstate]) == 1)
#     # Transition to a new true state
#     if (t < len)
#       trueState <- rmulti(1, size = 1, prob = probTrans[trueState, , t])
#   }
#   return(ans)
# })

#' @export
#' @rdname dDHMM
rDHMMo <- nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 init = double(1),
                 probObs = double(3),
                 probTrans = double(3),
                 len = double()) {
  returnType(double(1))
  ans <- numeric(len)

  trueInit <- 0

  r <- runif(1, 0, 1)
  j <- 1
  while (r > sum(init[1:j])) j <- j + 1
  trueState <- j

  for (i in 1:len) {
    # Detect based on the true state
    r <- runif(1, 0, 1)
    j <- 1
    while (r > sum(probObs[trueState, 1:j, i])) j <- j + 1
    ans[i] <- j

    # Transition to a new true state
    if (i != len) {
      r <- runif(1, 0, 1)
      j <- 1
      while (r > sum(probTrans[trueState, 1:j, i])) j <- j + 1
      trueState <- j
    }
  }
  return(ans)
})
# rDHMMo <- nimbleFunction(
#   run = function(n = integer(),    ## Observed capture (state) history
#                  init = double(1), ## probabilities of state at time 1
#                  probObs = double(3),
#                  probTrans = double(3),
#                  len = integer()) {
#     returnType(double(1))
#     ans <- numeric(len)
#
#     trueState <- rmulti(1, size = 1, prob = init)
#
#     for (t in 1:len) {
#       thisstate <- which(trueState == 1)
#       ans[t] <- which(rmulti(1, size = 1, prob = probObs[, thisstate, t]) == 1)
#       # Transition to a new true state
#       if (t < len)
#         trueState <- rmulti(1, size = 1, prob = probTrans[trueState, , t])
#     }
#     return(ans)
#   })
