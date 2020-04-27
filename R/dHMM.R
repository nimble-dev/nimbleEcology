#' Hidden Markov Model distribution for use in \code{nimble} models
#'
#' \code{dHMM} and \code{dHMMo} provide hidden Markov model
#' distributions that can be used directly from R or in \code{nimble}
#' models.
#'
#' @name dHMM
#'
#' @aliases dHMM dHMMo rHMM rHMMo
#'
#' @author Ben Goldstein, Perry de Valpine, and Daniel Turek
#'
#' @param x vector of observations, each one a positive integer
#'     corresponding to an observation state
#'     (one value of which could can correspond to "not observed", and
#'     another value of which can correspond to "dead" or
#'     "removed from system").
#' @param init vector of initial state probabilities. Must sum to 1
#' @param probObs time-independent matrix (\code{dHMM} and
#'     \code{rHMM}) or time-dependent array (\code{dHMMo} and
#'     \code{rHMMo}) of observation probabilities.
#'     First two dimensions of \code{probObs} are of size x (number of possible
#'     system states) x (number of possible observation classes). \code{dHMMo} and
#'     \code{rHMMo}
#'     expects an additional third dimension of size (number of observation times).
#' @param probTrans time-independent matrix of state transition
#'     probabilities.
#' @param len length of \code{x} (see below).
#' @param log TRUE or 1 to return log probability. FALSE or 0 to
#'     return probability.
#' @param n number of random draws, each returning a vector of length
#'     \code{len}. Currently only \code{n = 1} is supported, but the
#'     argument exists for standardization of "\code{r}" functions.
#'
#' @details
#' These nimbleFunctions provide distributions that can be used directly in R or
#' in \code{nimble} hierarchical models (via \code{\link[nimble]{nimbleCode}}
#' and \code{\link[nimble]{nimbleModel}}).
#'
#' The distribution has two forms, \code{dHMM} and
#' \code{dHMMo}. Define S as the number of latent state categories
#' (maximum possible value for elements of \code{x}), O as the number
#' of possible observation state categories, and T as the number of
#' observation times (length of \code{x}). In \code{dHMM},
#' \code{probObs} is a time-independent observation probability matrix
#' with dimension S x O.  In \code{dHMMo}, \code{probObs} is a
#' three-dimensional array of time-dependent observation probabilities
#' with dimension S x O x T.  The first index of \code{probObs}
#' indexes the true latent state.  The second index of \code{probObs}
#' indexes the observed state.  For example, in the time-dependent
#' case, \code{probObs[i, j, t]} is the probability at time \code{t} that
#' an individual in state \code{i} is observed in state \code{j}.
#'
#' \code{probTrans} has dimension S x S. \code{probTrans}[i, j] is the
#' tune-independent probability that an individual in state \code{i} at
#' time \code{t} transitions to state \code{j} time \code{t+1}.
#'
#' \code{initStates} has length S. \code{initStates[i]} is the
#' probability of being in state \code{i} at the first observation time.
#'
#' For more explanation, see
#' \href{../doc/Introduction_to_nimbleEcology.html}{package vignette} (or
#' \code{vignette("Introduction_to_nimbleEcology")}).
#'
#' Compared to writing \code{nimble} models with a discrete latent
#' state and a separate scalar datum for each observation time, use of
#' these distributions allows one to directly sum (marginalize) over
#' the discrete latent state and calculate the probability of all
#' observations for one individual (or other HMM unit) jointly.
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
#' \code{observedStates[i, 1:T] ~ dHMM(initStates[1:S],
#' observationProbs[1:S, 1:O],
#' transitionProbs[1:S, 1:S], T)}
#'
#' declares that the \code{observedStates[i, 1:T]} (observation
#' history for individual \code{i}, for example) vector follows a
#' hidden Markov model distribution with parameters as indicated,
#' assuming all the parameters have been declared elsewhere in the
#' model. As above, \code{S} is the number of system state categories,
#' \code{O} is the number of observation state categories, and
#' \code{T} is the number of observation occasions. This will invoke
#' (something like) the following call to \code{dHMM} when
#' \code{nimble} uses the model such as for MCMC:
#'
#' \code{dHMM(observedStates[1:T], initStates[1:S],
#' observationProbs[1:S, 1:O],
#' transitionProbs[1:S, 1:S], T, log = TRUE)}
#'
#' If an algorithm using a \code{nimble} model with this declaration
#' needs to generate a random draw for \code{observedStates[1:T]}, it
#' will make a similar invocation of \code{rHMM}, with \code{n = 1}.
#'
#' If the observation probabilities are time-dependent, one would use:
#'
#' \code{observedStates[1:T] ~ dHMMo(initStates[1:S], observationProbs[1:O, 1:S, 1:T],
#' transitionProbs[1:S, 1:S], T)}
#'
#' @return
#' For \code{dHMM} and \code{dHMMo}: the probability (or likelihood) or log
#' probability of observation vector \code{x}.
#'
#' For \code{rHMM} and \code{rHMMo}: a simulated detection history, \code{x}.
#'
#' @seealso For dynamic hidden Markov models with time-dependent transitions,
#' see \code{\link{dDHMM}} and \code{\link{dDHMMo}}.
#' For simple capture-recapture, see \code{\link{dCJS}}.
#'
#' @references D. Turek, P. de Valpine and C. J. Paciorek. 2016. Efficient Markov chain Monte
#' Carlo sampling for hierarchical hidden Markov models. Environmental and Ecological Statistics
#' 23:549–564. DOI 10.1007/s10651-016-0353-z
#'
#' @examples
#' \donttest{
#' # Set up constants and initial values for defining the model
#' len <- 5 # length of dataset
#' dat <- c(1,2,1,1,2) # A vector of observations
#' init <- c(0.4, 0.2, 0.4) # A vector of initial state probabilities
#' probObs <- t(array( # A matrix of observation probabilities
#'        c(1, 0,
#'          0, 1,
#'          0.2, 0.8), c(2, 3)))
#' probTrans <- t(array( # A matrix of transition probabilities
#'         c(0.6, 0.3, 0.1,
#'           0, 0.7, 0.3,
#'           0, 0, 1), c(3,3)))
#'
#' # Define code for a nimbleModel
#'  nc <- nimbleCode({
#'    x[1:5] ~ dHMM(init[1:3], probObs = probObs[1:3,1:2],
#'                  probTrans = probTrans[1:3, 1:3], len = 5)
#'
#'    for (i in 1:3) {
#'      for (j in 1:3) {
#'        probTrans[i,j] ~ dunif(0,1)
#'      }
#'
#'      probObs[i, 1] ~ dunif(0,1)
#'      probObs[i, 2] <- 1 - probObs[i, 1]
#'    }
#'  })
#'
#' # Build the model
#' HMM_model <- nimbleModel(nc,
#'                          data = list(x = dat),
#'                          inits = list(init = init,
#'                                       probObs = probObs,
#'                                       probTrans = probTrans))
#' # Calculate log probability of data from the model
#' HMM_model$calculate()
#' # Use the model for a variety of other purposes...
#' }

NULL

#' @export
#' @rdname dHMM
dHMM <- nimbleFunction(
  run = function(x = double(1), ## Observed capture (state) history
                 init = double(1),
                 probObs = double(2),
                 probTrans = double(2),
                 len = double(0, default = 0),## length of x (needed as a separate param for rDHMM)
                 log = integer(0, default = 0)) {
    if (length(x) != len) stop("Argument len must be length of x or 0.")
    if (dim(probObs)[1] != dim(probTrans)[1]) stop("Number of cols in probObs must equal number of cols in probTrans.")
    if (dim(probTrans)[1] != dim(probTrans)[2]) stop("probTrans must be a square matrix.")
    if (sum(init) != 1) stop("Initial probabilities must sum to 1.")

    pi <- init # State probabilities at time t=1
    logL <- 0
    nStates <- dim(probObs)[2]
    for (t in 1:len) {
      if (x[t] > nStates) stop("Invalid value of x[t] in dDHMM.")
      Zpi <- probObs[, x[t]] * pi # Vector of P(state) * P(observation class x[t] | state)
      sumZpi <- sum(Zpi)    # Total P(observed as class x[t])
      logL <- logL + log(sumZpi)  # Accumulate log probabilities through time
      if (t != len) pi <- (probTrans[,] %*% asCol(Zpi) / sumZpi)[ ,1] # State probabilities at t+1
    }
    returnType(double())
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @export
#' @rdname dHMM
dHMMo <- nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 init = double(1),##
                 probObs = double(3),
                 probTrans = double(2),
                 len = double(0, default = 0),## length of x (needed as a separate param for rDHMM)
                 log = integer(0, default = 0)) {
    if (length(x) != len) stop("Argument len must be length of x or 0.")
    if (dim(probObs)[1] != dim(probTrans)[1]) stop("Number of cols in Z must equal number of cols in T.")
    if (dim(probTrans)[1] != dim(probTrans)[2]) stop("T must be a square matrix.")
    if (dim(probObs)[3] != len) {
      if (dim(probObs)[3] == 1) stop("Time dimension of Z must match length of data. Did you mean dHMM?")
      stop("Length of time dimension of Z must match length of data.")
    }
    if (sum(init) != 1) stop("Initial probabilities must sum to 1.")

    pi <- init # State probabilities at time t=1
    logL <- 0
    nStates <- dim(probObs)[2]
    for (t in 1:len) {
      if (x[t] > nStates) stop("Invalid value of x[t] in dDHMM.")
      Zpi <- probObs[,x[t],t] * pi # Vector of P(state) * P(observation class x[t] | state)
      sumZpi <- sum(Zpi)    # Total P(observed as class x[t])
      logL <- logL + log(sumZpi)  # Accumulate log probabilities through time
      if (t != len) pi <- (probTrans[,] %*% asCol(Zpi) / sumZpi)[ ,1] # State probabilities at t+1
    }
    returnType(double())
    if (log) return(logL)
    return(exp(logL))
  }
)

#' @export
#' @rdname dHMM
rHMM <- nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 init = double(1),
                 probObs = double(2),
                 probTrans = double(2),
                 len = double(0, default = 0)) {
  returnType(double(1))
  ans <- numeric(len)

  probInit <- init
  trueInit <- 0

  r <- runif(1, 0, 1)
  j <- 1
  while (r > sum(probInit[1:j])) j <- j + 1
  trueInit <- j

  trueState <- trueInit
  ### QUESTION: Is the "init" probability for the state at time t1 or t0? I'm assuming t0
  for (i in 1:len) {
    # Transition to a new true state
    r <- runif(1, 0, 1)
    j <- 1
    while (r > sum(probTrans[trueState, 1:j])) j <- j + 1
    trueState <- j

    # Detect based on the true state
    r <- runif(1, 0, 1)
    j <- 1
    while (r > sum(probObs[trueState, 1:j])) j <- j + 1
    ans[i] <- j

  }

  return(ans)
})

#' @export
#' @rdname dHMM
rHMMo <- nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 init = double(1),
                 probObs = double(3),
                 probTrans = double(2),
                 len = double(0, default = 0)) {
  returnType(double(1))
  ans <- numeric(len)

  probInit <- init
  trueInit <- 0

  r <- runif(1, 0, 1)
  j <- 1
  while (r > sum(probInit[1:j])) j <- j + 1
  trueInit <- j

  trueState <- trueInit
  for (i in 1:len) {
    # Transition to a new true state
    r <- runif(1, 0, 1)
    j <- 1
    while (r > sum(probTrans[trueState, 1:j])) j <- j + 1
    trueState <- j

    # Detect based on the true state
    r <- runif(1, 0, 1)
    j <- 1
    while (r > sum(probObs[trueState, 1:j, i])) j <- j + 1
    ans[i] <- j

  }

  return(ans)
})



