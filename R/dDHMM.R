#' Dynamic Hidden Markov Model distribution for use in NIMBLE models
#'
#' \code{dDHMM} and \code{dDHMMo} provide Dynamic hidden Markov model distributions for NIMBLE models.
#' "Dynamic" here means that the matrix of state transition probabilities in indexed by time.  The
#' \code{dDHMMo} alias is used when observation probabilities are indexed by time.
#'
#' Compared to writing NIMBLE models with discrete latent states, use of these DHMM distributions allows
#' one to directly integrate over such discrete latent states and hence leave them out of the NIMBLE
#' model code.
#'
#' @aliases dDHMM dDHMMo rDHMM rDHMMo
#'
#' @export
#'
#' @name dDHMM
#'
#' @param x vector of observation classes, one of which could be defined as "not observed".
#' @param init vector of initial state probabilities
#' @param Z time-independent matrix of observation probabilities.
#' First two dimensions of \code{Z} are of size (number of possible observation classes) x
#'  (number of possible system states).  In \code{dDHMMo}, the third dimension of \code{Z} is of
#'  size (number of observation times).
#' @param T time-dependent matrix of system state-transition probabilities.
#' Dimension of \code{T} is (number of possible system states) x  (number of possible system states)
#' x (number of observation times).
#' @param len length of observations (needed for rDHMM)
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return probability.
#' @param n length of random sequence
#'
#' @author Perry de Valpine, Daniel Turek, and Ben Goldstein
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
#' \code{observedStates[1:T] ~ dDHMM(initStates[1:S], observationProbs[1:O, 1:S],
#' transitionProbs[1:S, 1:S, 1:T], T)}
#'
#' declares that the \code{observedStates[1:T]} vector follows a dynamic hidden Markov model distribution
#' with parameters as indicated, assuming all the parameters have been declared elsewhere in the model.  In
#' this case, \code{S} is the number of system states, \code{O} is the number of observation classes, and
#' \code{T} is the number of observation occasions.
#'
#' If the observation probabilities are time-dependent, one would use:
#'
#' \code{observedStates[1:T] ~ dDHMMo(initStates[1:S], observationProbs[1:O, 1:S, 1:T],
#' transitionProbs[1:S, 1:S, 1:T], T)}
#'
#' @seealso For hidden Markov models with time-independent transitions, see \link{dHMM} and \link{dHMMo}.
#' For simple capture-recapture, see \link{dCJS}.

NULL

#' @export
#' @rdname dDHMM
dDHMM <- nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 init = double(1),
                 Z = double(2),
                 T = double(3),
                 len = double(),## length of x (needed as a separate param for rDHMM)
                 log = integer(0, default = 0)) {
    if (length(init) != dim(Z)[2]) stop("Length of init does not match ncol of Z in dDHMM.")
    if (length(init) != dim(T)[1]) stop("Length of init does not match dim(T)[1] in dDHMM.")
    if (length(init) != dim(T)[2]) stop("Length of init does not match dim(T)[2] in dDHMM.")
    if (length(x) != len) stop("Length of x does not match len in dDHMM.")
    if (len-1 != dim(T)[3]) stop("Length of x does not match len in dDHMM.")

    pi <- init # State probabilities at time t=1
    logL <- 0
    nObsClasses <- dim(Z)[1]
    lengthX <- length(x)
    for (t in 1:lengthX) {
      if (x[t] > nObsClasses) stop("Invalid value of x[t] in dDHMM.")
      Zpi <- Z[x[t], ] * pi # Vector of P(state) * P(observation class x[t] | state)
      sumZpi <- sum(Zpi)    # Total P(observed as class x[t])
      logL <- logL + log(sumZpi)  # Accumulate log probabilities through time
      if (t != lengthX)   pi <- (T[,,t] %*% asCol(Zpi) / sumZpi)[ ,1] # State probabilities at t+1
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
                 Z = double(3),
                 T = double(3),
                 len = double(),## length of x (needed as a separate param for rDHMM)
                 log = integer(0, default = 0)) {
    if (length(init) != dim(Z)[2]) stop("Length of init does not match ncol of Z in dDHMMo.")
    if (length(init) != dim(T)[1]) stop("Length of init does not match dim(T)[1] in dDHMMo.")
    if (length(init) != dim(T)[2]) stop("Length of init does not match dim(T)[2] in dDHMMo.")
    if (length(x) != len) stop("Length of x does not match len in dDHMM.")
    if (len != dim(T)[3]) stop("dim(T)[3] does not match len in dDHMMo.")
    if (len != dim(Z)[3]) stop("dim(Z)[3] does not match len in dDHMMo.")

    pi <- init # State probabilities at time t=1
    logL <- 0
    nObsClasses <- dim(Z)[1]
    lengthX <- length(x)
    for (t in 1:lengthX) {
      if (x[t] > nObsClasses) stop("Invalid value of x[t] in dDHMM.")
      Zpi <- Z[x[t], , t] * pi # Vector of P(state) * P(observation class x[t] | state)
      sumZpi <- sum(Zpi)    # Total P(observed as class x[t])
      logL <- logL + log(sumZpi)  # Accumulate log probabilities through time
      if (t != lengthX)   pi <- (T[,,t] %*% asCol(Zpi) / sumZpi)[ ,1] # State probabilities at t+1
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
                 Z = double(2),
                 T = double(3),
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
    while (r > sum(T[trueState, 1:j, i])) j <- j + 1
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
#' @rdname dDHMM
rDHMMo <- nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 init = double(1),
                 Z = double(3),
                 T = double(3),
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
    while (r > sum(T[trueState, 1:j, i])) j <- j + 1
    trueState <- j

    # Detect based on the true state
    r <- runif(1, 0, 1)
    j <- 1
    while (r > sum(Z[1:j, trueState, i])) j <- j + 1
    ans[i] <- j
  }

  return(ans)
})
