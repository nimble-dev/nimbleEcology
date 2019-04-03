#' Hidden Markov Model distribution for use in NIMBLE models
#'
#' \code{dHMM} and \code{dHMMo} provide hidden Markov model distributions for NIMBLE models.
#' "Dynamic" here means that the matrix of state transition probabilities is indexed by time.  The
#' \code{dHMMo} version additionally allows observation probabilities to be indexed by time.
#' Compared to writing NIMBLE models with discrete latent states, use of these DHMM distributions allows
#' one to directly integrate over such discrete latent states and hence leave them out of the NIMBLE
#' model code.
#'
#' @aliases dHMM dHMMo rHMM rHMMo
#'
#' @name dHMM
#'
#' @export
#'
#' @param x vector of observation classes, one of which could be defined as "not observed".
#' @param init vector of initial state probabilities
#' @param Z time-independent matrix of observation probabilities.
#' Dimension of \code{Z} is (number of possible observation classes) x (number of possible system states)
#' @param T time-independent matrix of system state-transition probabilities.
#' Dimension of \code{T} is (number of possible system states) x (number of possible system states)
#' @param len length of observations (needed for rDHMM)
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return probability.
#'
#' @author Ben Goldstein, Perry de Valpine, and Daniel Turek
#'
#' @references D. Turek, P. de Valpine and C. J. Paciorek. 2016. Efficient Markov chain Monte
#' Carlo sampling for hierarchical hidden Markov models. Environmental and Ecological Statistics
#' 23:549–564. DOI 10.1007/s10651-016-0353-z
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
#' \code{observedStates[1:T] ~ dHMM(initStates[1:S], observationProbs[1:O, 1:S],
#' transitionProbs[1:S, 1:S], T)}
#'
#' declares that the \code{observedStates[1:T]} vector follows a hidden Markov model distribution
#' with parameters as indicated, assuming all the parameters have been declared elsewhere in the model.  In
#' this case, \code{S} is the number of system states and \code{O} is the number of observation classes, and
#' \code{T} is the number of observation occasions.
#'
#' If the observation probabilities are time-dependent, one would use:
#'
#' \code{observedStates[1:T] ~ dHMMo(initStates[1:S], observationProbs[1:O, 1:S, 1:T],
#' transitionProbs[1:S, 1:S], T)}
#'
#' @seealso For dynamic hidden Markov models with time-dependent transitions, see \link{dDHMM} and \link{dDHMMo}.
#' For simple capture-recapture, see \link{dCJS}.

#' @export
#' @rdname dHMM
dHMM <- nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 init = double(1),##
                 Z = double(2),
                 T = double(2),
                 len = double(0, default = 0),## length of x (needed as a separate param for rDHMM)
                 log = integer(0, default = 0)) {
    if (length(x) != len) stop("Argument len must be length of x or 0.")
    if (dim(Z)[2] != dim(T)[1]) stop("Number of cols in Z must equal number of cols in T.")
    if (dim(T)[1] != dim(T)[2]) stop("T must be a square matrix.")
    if (sum(init) != 1) stop("Initial probabilities must sum to 1.")

    pi <- init # State probabilities at time t=1
    logL <- 0
    nStates <- dim(Z)[1]
    for (t in 1:len) {
      if (x[t] > nStates) stop("Invalid value of x[t] in dDHMM.")
      Zpi <- Z[x[t], ] * pi # Vector of P(state) * P(observation class x[t] | state)
      sumZpi <- sum(Zpi)    # Total P(observed as class x[t])
      logL <- logL + log(sumZpi)  # Accumulate log probabilities through time
      if (t != len) pi <- (T[,] %*% asCol(Zpi) / sumZpi)[ ,1] # State probabilities at t+1
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
                 Z = double(3),
                 T = double(2),
                 len = double(0, default = 0),## length of x (needed as a separate param for rDHMM)
                 log = integer(0, default = 0)) {
    if (length(x) != len) stop("Argument len must be length of x or 0.")
    if (dim(Z)[2] != dim(T)[1]) stop("Number of cols in Z must equal number of cols in T.")
    if (dim(T)[1] != dim(T)[2]) stop("T must be a square matrix.")
    if (dim(Z)[3] != len) {
      if (dim(Z)[3] == 1) stop("Time dimension of Z must match length of data. Did you mean dHMM?")
      stop("Length of time dimension of Z must match length of data.")
    }
    if (sum(init) != 1) stop("Initial probabilities must sum to 1.")

    pi <- init # State probabilities at time t=1
    logL <- 0
    nStates <- dim(Z)[1]
    for (t in 1:len) {
      if (x[t] > nStates) stop("Invalid value of x[t] in dDHMM.")
      Zpi <- Z[x[t],,t] * pi # Vector of P(state) * P(observation class x[t] | state)
      sumZpi <- sum(Zpi)    # Total P(observed as class x[t])
      logL <- logL + log(sumZpi)  # Accumulate log probabilities through time
      if (t != len) pi <- (T[,] %*% asCol(Zpi) / sumZpi)[ ,1] # State probabilities at t+1
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
                 Z = double(2),
                 T = double(2),
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
    while (r > sum(T[trueState, 1:j])) j <- j + 1
    trueState <- j

    # Detect based on the true state
    r <- runif(1, 0, 1)
    j <- 1
    while (r > sum(Z[1:j, trueState])) j <- j + 1
    ans[i] <- j

  }

  return(ans)
})

#' @export
#' @rdname dHMM
rHMMo <- nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 init = double(1),
                 Z = double(3),
                 T = double(2),
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
    while (r > sum(T[trueState, 1:j])) j <- j + 1
    trueState <- j

    # Detect based on the true state
    r <- runif(1, 0, 1)
    j <- 1
    while (r > sum(Z[1:j, trueState, i])) j <- j + 1
    ans[i] <- j

  }

  return(ans)
})


registerDistributions(list(
  dHMM = list(
    BUGSdist = "dHMM(init, Z, T, len)",
    Rdist = "dHMM(init, Z, T, len = 0)",
    discrete = TRUE,
    types = c('value = double(1)', 'init = double(1)', 'Z = double(2)', 'T = double(2)', 'len = double(0)'),
    pqAvail = FALSE))
  )

registerDistributions(list(
  dHMMo = list(
    BUGSdist = "dHMMo(init, Z, T, len)",
    Rdist = "dHMMo(init, Z, T, len = 0)",
    discrete = TRUE,
    types = c('value = double(1)', 'init = double(1)', 'Z = double(3)', 'T = double(2)', 'len = double(0)'),
    pqAvail = FALSE))
  )

