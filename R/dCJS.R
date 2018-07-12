#' Cormack-Jolly-Seber distribution for use in NIMBLE models
#'
#' \code{dCJSss} and \code{dCJSvv} provide a basic distribution for capture history vectors based
#'  on survival and capture probabilities.  The different names are for scalar ("s", time-independent)
#'  versus vector ("v", time-dependent) survival and capture probabilities.
#'
#' @aliases dCJSvv rCJSss rCJSvv dCJS
#'
#' @export
#'
#' @param x capture-history vector of 0s (not captured) and 1s (captured), assumed to begin with a 1.
#' @param probSurvive survival probability, either a scalar (for dCJSss) or a vector (for dCJSvv).
#' @param probCapture capture probability, either a scalar (for dCJSss) or a vector (for dCJSvv).
#' @param l length of capture history (needed for rCJSxx).
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return probability.
#'
#' @author Perry de Valpine and Daniel Turek
#'
#' @details These nimbleFunctions provide distributions that can be used in code (via \code{nimbleCode}) for \link{nimbleModel}.
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
#' For example,
#'
#' \code{captures[1:T] ~ dCSJss(survive, capture, T)}
#'
#' declares a vector node, \code{captures[1:T]}, that follows a capture-recapture distribution
#' with scalar survival probability \code{survive} and scalar capture probability \code{capture}
#' (assuming \code{survive} and \code{capture} are defined elsewhere in the model).
#' If time-dependent survival and capture probabilities are needed, use
#'
#' \code{captures[1:T] ~ dCSJss(survive[1:T], capture[1:T], T)}.
#'
#' In fact, the \code{len} argument (\code{T} in these examples) will be ignored during model
#' probability calculations.  It will be used only for simulations.
#'
#' \emph{Ignored elements of \code{survive} and \code{capture} when they are vectors:} Note that the last element of \code{survive},
#' if present, will be ignored.  Also the first element of \code{capture} will be ignored, since the first
#' element of the capture history is assumed to be 1, indicated first capture.
#'
#' @references D. Turek, P. de Valpine and C. J. Paciorek. 2016. Efficient Markov chain Monte
#' Carlo sampling for hierarchical hidden Markov models. Environmental and Ecological Statistics
#' 23:549–564. DOI 10.1007/s10651-016-0353-z
#'
#' @seealso For multi-state or multi-event capture-recapture models, see \link{dHMM} or \link{dDHMM}.
dCJSss <- nimbleFunction(
  run = function(x = double(1),    ## standard name for the "data"
                 probSurvive = double(),
                 probCapture = double(),
                 len = double(),
                 log = integer(0, default = 0) ## required log argument
  ) {
    ## Note the calculations used here are actually in hidden Markov model form.
    probAliveGivenHistory <- 1
    ## logProbData will be the final answer
    logProbData <- 0
    if(len == 1) {  ## l==1 should not occur, but just in case:
      return(logProbData)
    }
    for(t in 2:len) {
      ## probAlive is P(Alive(t) | x(1)...x(t-1))
      ## probAliveGivenHistory is (Alive(t-1) | x(1)...x(t-1))
      probAlive <- probAliveGivenHistory * probSurvive
      if(x[t] == 1) {
        ## ProbThisObs = P(x(t) | x(1)...x(t-1))
        probThisObs <- probAlive * probCapture
        probAliveGivenHistory <- 1
      } else {
        probAliveNotSeen <- probAlive * (1-probCapture)
        probThisObs <- probAliveNotSeen + (1-probAlive)
        probAliveGivenHistory <- probAliveNotSeen / probThisObs
      }
      logProbData <- logProbData + log(probThisObs)
    }
    if(log) return(logProbData)
    return(exp(logProbData))
    returnType(double())
  }
)

dCJSvv <- nimbleFunction(
  # The first element of x is ignored.
  # It is assumed to represent first capture and hence to be 1.
  # probSurvive[t] represents survival from t to t+1, so the last element of probSurvive is ignored.
  # probCapture[t] represents capture probability at time t, so the first element of probCapture is ignored.
  run = function(x = double(1),    ## standard name for the "data"
                 probSurvive = double(1),
                 probCapture = double(1),
                 len = double(),
                 log = integer(0, default = 0) ## required log argument
  ) {
    ## Note the calculations used here are actually in hidden Markov model form.
    probAliveGivenHistory <- 1
    ## logProbData will be the final answer
    logProbData <- 0
    if(len == 1) {  ## l==1 should not occur, but just in case:
      return(logProbData)
    }
    for(t in 2:len) {
      ## probAlive is P(Alive(t) | x(1)...x(t-1))
      ## probAliveGivenHistory is (Alive(t-1) | x(1)...x(t-1))
      probAlive <- probAliveGivenHistory * probSurvive[1]
      if(x[t] == 1) {
        ## ProbThisObs = P(x(t) | x(1)...x(t-1))
        probThisObs <- probAlive * probCapture[i]
        probAliveGivenHistory <- 1
      } else {
        probAliveNotSeen <- probAlive * (1-probCapture[i])
        probThisObs <- probAliveNotSeen + (1-probAlive)
        probAliveGivenHistory <- probAliveNotSeen / probThisObs
      }
      logProbData <- logProbData + log(probThisObs)
    }
    if(log) return(logProbData)
    return(exp(logProbData))
    returnType(double())
  }
)

rCJSvv <- nimbleFunction(
  run = function(n = integer(),
                 probSurvive = double(1),
                 probCapture = double(1),
                 len = double()) {
    ans <- numeric(len)
    ans[1] <- 1
    alive <- 1
    if(l == 1) return(ans)
    for(i in 2:len) {
      if(alive)
        alive <- rbinom(1, size = 1, prob = probSurvive[i])
      if(alive) {
        ans[i] <- rbinom(1, size = 1, prob = probCapture[i])
      } else {
        ans[i] <- 0
      }
    }
    return(ans)
    returnType(double(1))
  }
)


rCJSss <- nimbleFunction(
  run = function(n = integer(),
                 probSurvive = double(),
                 probCapture = double(),
                 len = double()) {
    ans <- numeric(len)
    ans[1] <- 1
    alive <- 1
    if(l == 1) return(ans)
    for(i in 2:len) {
      if(alive)
        alive <- rbinom(1, size = 1, prob = probSurvive)
      if(alive) {
        ans[i] <- rbinom(1, size = 1, prob = probCapture)
      } else {
        ans[i] <- 0
      }
    }
    return(ans)
    returnType(double(1))
  }
)
