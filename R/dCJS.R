#' Cormack-Jolly-Seber distribution for use in NIMBLE models
#'
#' \code{dCJS_**} functions provide a basic distribution for capture history vectors based
#'  on survival and capture probabilities. The different aliases are for scalar ("s", time-independent)
#'  versus vector ("v", time-dependent) survival and capture probabilities, in that order.
#'
#' @aliases dCJS_ss dCJS_sv dCJS_vs dCJS_vv rCJS_ss rCJS_sv rCJS_vs rCJS_vv
#'
#' @name dCJS
#'
#' @param x capture-history vector of 0s (not captured) and 1s (captured). Do not include the initial capture, which is assumed.
#' @param probSurvive survival probability, either a scalar (for dCJS_s*) or a vector (for dCJS_v*).
#' @param probCapture capture probability, either a scalar (for dCJS_*s) or a vector (for dCJS_*v).
#' @param len length of capture history (needed for rCJSxx).
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return probability.
#' @param n length of random sequence
#'
#' @author Ben Goldstein, Perry de Valpine, and Daniel Turek
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
#' \code{captures[1:T] ~ dCSJ_ss(survive, capture, T)}
#'
#' declares a vector node, \code{captures[1:T]}, that follows a capture-recapture distribution
#' with scalar survival probability \code{survive} and scalar capture probability \code{capture}
#' (assuming \code{survive} and \code{capture} are defined elsewhere in the model).
#' If time-dependent survival and capture probabilities are needed, use
#'
#' \code{captures[1:T] ~ dCSJ_vv(survive[1:T], capture[1:T], T)}.
#'
#' In fact, the \code{len} argument (\code{T} in these examples) will be ignored during model
#' probability calculations.  It will be used only for simulations.
#'
#' It is important to use the time indexing correctly for survival.  \code{probSurvive[t]} is the survival
#' probabilty from time \code{t-1} to time \code{t}.  Time indexing for detection is more obvious:
#' \code{probDetect[t]} is the detection probability at time \code{t}.
#'
#' @references D. Turek, P. de Valpine and C. J. Paciorek. 2016. Efficient Markov chain Monte
#' Carlo sampling for hierarchical hidden Markov models. Environmental and Ecological Statistics
#' 23:549–564. DOI 10.1007/s10651-016-0353-z
#'
#' @seealso For multi-state or multi-event capture-recapture models, see \link{dHMM} or \link{dDHMM}.
#' @import nimble
#' @importFrom stats rbinom runif dbinom
#'
#' @examples
#' \dontrun{
#' # Set up constants and initial values for defining the model
#' dat <- c(1,1,0,0) # A vector of observations
#' probSurvive <- 0.6
#' probCapture <- 0.4
#'
#'
#' # Define code for a nimbleModel
#' nc <- nimbleCode({
#'   x[1:4] ~ dCJSss(probSurvive, probCapture, len = 4)
#'   probSurvive ~ dunif(0,1)
#'   probCapture ~ dunif(0,1)
#' })
#'
#' # Build the model, providing data and initial values
#' CJS_model <- nimbleModel(nc, data = list(x = dat),
#'                          inits = list(probSurvive = probSurvive,
#'                                       probCapture = probCapture))
#'
#' # Calculate log probability of data from the model
#' CJS_model$calculate()
#' # Use the model for a variety of other purposes...
#' }


NULL

#' @rdname dCJS
#' @export
dCJS_ss <- nimbleFunction(
  run = function(x = double(1),    ## standard name for the "data"
                 probSurvive = double(),
                 probCapture = double(),
                 len = double(0, default = 0),
                 log = integer(0, default = 0) ## required log argument
  ) {

    if (len != 0) {
      if (len != length(x)) stop("Argument len must match length of data, or be 0.")
    }

    ## Note the calculations used here are actually in hidden Markov model form.
    probAliveGivenHistory <- 1
    ## logProbData will be the final answer
    logProbData <- 0
    lenX <- length(x)
    if (lenX == 0) {  ## lenX < 1 should not occur, but just in case:
      return(0)
    }
    for (t in 1:lenX) {
      ## probAlive is P(Alive(t) | x(1)...x(t-1))
      ## probAliveGivenHistory is (Alive(t-1) | x(1)...x(t-1))
      probAlive <- probAliveGivenHistory * probSurvive
      if (x[t] == 1) {
        ## ProbThisObs = P(x(t) | x(1)...x(t-1))
        probThisObs <- probAlive * probCapture
        probAliveGivenHistory <- 1
      } else {
        probAliveNotSeen <- probAlive * (1 - probCapture)
        probThisObs <- probAliveNotSeen + (1 - probAlive)
        probAliveGivenHistory <- probAliveNotSeen / probThisObs
      }
      logProbData <- logProbData + log(probThisObs)
    }
    if (log) return(logProbData)
    return(exp(logProbData))
    returnType(double(0))
  }
)

#' @rdname dCJS
#' @export
dCJS_sv <- nimbleFunction(
  run = function(x = double(1),    ## standard name for the "data"
                 probSurvive = double(),
                 probCapture = double(1),
                 len = double(0, default = 0),
                 log = integer(0, default = 0) ## required log argument
  ) {
    if (len != 0) {
      if (len != length(x)) stop("Argument len must match length of data, or be 0.")
    }
    if (length(x) != length(probCapture)) stop("Length of probCapture does not match length of data.")

    ## Note the calculations used here are actually in hidden Markov model form.
    probAliveGivenHistory <- 1
    ## logProbData will be the final answer
    logProbData <- 0
    lenX <- length(x)
    if (lenX == 0) {  ## l < 1 should not occur, but just in case:
      return(0)
    }
    for (t in 1:lenX) {
      ## probAlive is P(Alive(t) | x(1)...x(t-1))
      ## probAliveGivenHistory is (Alive(t-1) | x(1)...x(t-1))
      probAlive <- probAliveGivenHistory * probSurvive
      if (x[t] == 1) {
        ## ProbThisObs = P(x(t) | x(1)...x(t-1))
        probThisObs <- probAlive * probCapture[t]
        probAliveGivenHistory <- 1
      } else {
        probAliveNotSeen <- probAlive * (1 - probCapture[t])
        probThisObs <- probAliveNotSeen + (1 - probAlive)
        probAliveGivenHistory <- probAliveNotSeen / probThisObs
      }
      logProbData <- logProbData + log(probThisObs)
    }
    if (log) return(logProbData)
    return(exp(logProbData))
    returnType(double())
  }
)


#' @rdname dCJS
#' @export
dCJS_vs <- nimbleFunction(
  run = function(x = double(1),    ## standard name for the "data"
                 probSurvive = double(1),
                 probCapture = double(),
                 len = double(0, default = 0),
                 log = integer(0, default = 0) ## required log argument
  ) {
    if (len != 0) {
      if (len != length(x)) stop("Argument len must match length of data, or be 0.")
    }
    if (length(x) != length(probSurvive)) stop("Length of probSurvive does not match length of data.")

    ## Note the calculations used here are actually in hidden Markov model form.
    probAliveGivenHistory <- 1
    ## logProbData will be the final answer
    logProbData <- 0
    lenX <- length(x)
    if (lenX == 0) {  ## l < 1 should not occur, but just in case:
      return(0)
    }
    for (t in 1:lenX) {
      ## probAlive is P(Alive(t) | x(1)...x(t-1))
      ## probAliveGivenHistory is (Alive(t-1) | x(1)...x(t-1))
      probAlive <- probAliveGivenHistory * probSurvive[t]
      if (x[t] == 1) {
        ## ProbThisObs = P(x(t) | x(1)...x(t-1))
        probThisObs <- probAlive * probCapture
        probAliveGivenHistory <- 1
      } else {
        probAliveNotSeen <- probAlive * (1 - probCapture)
        probThisObs <- probAliveNotSeen + (1 - probAlive)
        probAliveGivenHistory <- probAliveNotSeen / probThisObs
      }
      logProbData <- logProbData + log(probThisObs)
    }
    if (log) return(logProbData)
    return(exp(logProbData))
    returnType(double())
  }
)


#' @rdname dCJS
#' @export
dCJS_vv <- nimbleFunction(
  # It is assumed that the individual has already been captured.
  # Therefore, the first entry in x represents the first possible recapture event.
  # probSurvive[t] represents survival from t-1 to t.
  # probCapture[t] represents capture probability at time t.
  run = function(x = double(1),    ## standard name for the "data"
                 probSurvive = double(1),
                 probCapture = double(1),
                 len = double(0, default = 0),
                 log = integer(0, default = 0) ## required log argument
  ) {
    if (len != 0) {
      if (len != length(x)) stop("Argument len must match length of data, or be 0.")
    }
    if (length(x) != length(probSurvive)) stop("Length of probSurvive does not match length of data.")
    if (length(x) != length(probCapture)) stop("Length of probCapture does not match length of data.")
    ## Note the calculations used here are actually in hidden Markov model form.
    probAliveGivenHistory <- 1
    ## logProbData will be the final answer
    logProbData <- 0
    if (len == 0) {  ## l<1 should not occur, but just in case:
      len <- length(x)
    }
    for (t in 1:len) {
      ## probAlive is P(Alive(t) | x(1)...x(t-1))
      ## probAliveGivenHistory is (Alive(t-1) | x(1)...x(t-1))
      probAlive <- probAliveGivenHistory * probSurvive[t]
      if (x[t] == 1) {
        ## ProbThisObs = P(x(t) | x(1)...x(t-1))
        probThisObs <- probAlive * probCapture[t]
        probAliveGivenHistory <- 1
      } else {
        probAliveNotSeen <- probAlive * (1 - probCapture[t])
        probThisObs <- probAliveNotSeen + (1 - probAlive)
        probAliveGivenHistory <- probAliveNotSeen / probThisObs
      }
      logProbData <- logProbData + log(probThisObs)
    }
    if (log) {
      return(logProbData)
    }
    return(exp(logProbData))
    returnType(double())
  }
)

#' @rdname dCJS
#' @export
rCJS_ss <- nimbleFunction(
  run = function(n = integer(),
                 probSurvive = double(),
                 probCapture = double(),
                 len = double(0, default = 0)) {
    if (n != 1) stop("rCJS only works for n = 1")
    if(len < 0)
      stop("len must be non-negative.")
    ans <- numeric(length = len, init = FALSE)
    alive <- 1
    if (len <= 0) return(ans)
    for (i in 1:len) {
      if (alive)
        alive <- rbinom(1, size = 1, prob = probSurvive)
      if (alive) {
        ans[i] <- rbinom(1, size = 1, prob = probCapture)
      } else {
        ans[i] <- 0
      }
    }
    return(ans)
    returnType(double(1))
  }
)

#' @rdname dCJS
#' @export
rCJS_sv <- nimbleFunction(
  run = function(n = integer(),
                 probSurvive = double(0),
                 probCapture = double(1),
                 len = double(0, default = 0)) {
    if (n != 1) stop("rCJS only works for n = 1")
    if(len < 0)
      stop("len must be non-negative.")
    if(length(probCapture) != len)
      stop("Length of probCapture is not the same as len.")
    ans <- numeric(length = len, init = FALSE)
    alive <- 1
    if (len <= 0) return(ans)
    for (i in 1:len) {
      if (alive)
        alive <- rbinom(1, size = 1, prob = probSurvive)
      if (alive) {
        ans[i] <- rbinom(1, size = 1, prob = probCapture[i])
      } else {
        ans[i] <- 0
      }
    }
    return(ans)
    returnType(double(1))
  }
)

#' @export
#' @rdname dCJS
rCJS_vs <- nimbleFunction(
  run = function(n = integer(),
                 probSurvive = double(1),
                 probCapture = double(0),
                 len = double(0, default = 0)) {
    if (n != 1) stop("rCJS only works for n = 1")
    if(len < 0)
      stop("len must be non-negative.")
    if(length(probSurvive) != len)
      stop("Length of probSurvive is not the same as len.")
    ans <- numeric(length = len, init = FALSE)
    alive <- 1
    if (len <= 0) return(ans)
    for (i in 1:len) {
      if (alive)
        alive <- rbinom(1, size = 1, prob = probSurvive[i])
      if (alive) {
        ans[i] <- rbinom(1, size = 1, prob = probCapture)
      } else {
        ans[i] <- 0
      }
    }
    return(ans)
    returnType(double(1))
  }
)

#' @rdname dCJS
#' @export
rCJS_vv <- nimbleFunction(
  run = function(n = integer(),
                 probSurvive = double(1),
                 probCapture = double(1),
                 len = double(0, default = 0)) {
    if (n != 1) stop("rCJS only works for n = 1")
    if(len < 0)
      stop("len must be non-negative.")
    if(length(probSurvive) != len)
      stop("Length of probSurvive is not the same as len.")
    if(length(probCapture) != len)
      stop("Length of probCapture is not the same as len.")
    ans <- numeric(lenth = len, init = FALSE)
    alive <- 1
    if (len <= 0) return(ans)
    for (i in 1:len) {
      if (alive)
        alive <- rbinom(1, size = 1, prob = probSurvive[i])
      if (alive) {
        ans[i] <- rbinom(1, size = 1, prob = probCapture[i])
      } else {
        ans[i] <- 0
      }
    }
    return(ans)
    returnType(double(1))
  }
)




# Register the distributions explicitly for two reasons:
# 1. Avoid message to user about automatic registrations upon first use in a nimbleModel
# 2. Establish default len = 0 via reparameterization mechanism.
registerDistributions(list(
  dCJS_ss = list(
    BUGSdist = "dCJS_ss(probSurvive, probCapture, len)",
    Rdist = "dCJS_ss(probSurvive, probCapture, len = 0)",
    discrete = TRUE,
    types = c('value = double(1)', 'probSurvive = double(0)', 'probCapture = double(0)', 'len = double(0)'),
    pqAvail = FALSE))
  )

registerDistributions(list(
  dCJS_sv = list(
    BUGSdist = "dCJS_sv(probSurvive, probCapture, len)",
    Rdist = "dCJS_sv(probSurvive, probCapture, len = 0)",
    discrete = TRUE,
    types = c('value = double(1)', 'probSurvive = double(0)', 'probCapture = double(1)', 'len = double(0)'),
    pqAvail = FALSE))
  )

registerDistributions(list(
  dCJS_vs = list(
    BUGSdist = "dCJS_vs(probSurvive, probCapture, len)",
    Rdist = "dCJS_vs(probSurvive, probCapture, len = 0)",
    discrete = TRUE,
    types = c('value = double(1)', 'probSurvive = double(1)', 'probCapture = double(0)', 'len = double(0)'),
    pqAvail = FALSE))
  )

registerDistributions(list(
  dCJS_vv = list(
    BUGSdist = "dCJS_vv(probSurvive, probCapture, len)",
    Rdist = "dCJS_vv(probSurvive, probCapture, len = 0)",
    discrete = TRUE,
    types = c('value = double(1)', 'probSurvive = double(1)', 'probCapture = double(1)', 'len = double(0)'),
    pqAvail = FALSE))
)
