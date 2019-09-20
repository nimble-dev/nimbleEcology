#' Cormack-Jolly-Seber distribution for use in \code{nimble} models
#'
#' \code{dCJS_**} \code{rCJS_**} provide Cormack-Jolly-Seber capture-recapture
#' distributions that can be used directly from R or in \code{nimble}
#' models.
#'
#' @aliases dCJS_ss dCJS_sv dCJS_vs dCJS_vv rCJS_ss rCJS_sv rCJS_vs rCJS_vv
#'
#' @name dCJS
#'
#' @param x capture-history vector of 0s (not captured) and 1s (captured). Do not include the initial capture, which is assumed to have occurred prior to \code{x[1]}.
#' @param probSurvive survival probability, either a scalar (for dCJS_s*) or a vector (for dCJS_v*).
#' @param probCapture capture probability, either a scalar (for dCJS_*s) or a vector (for dCJS_*v).
#' @param len length of capture history (see below).
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return probability.
#' @param n number of random draws, each returning a vector of length
#'     \code{len}. Currently only \code{n = 1} is supported, but the
#'     argument exists for standardization of "\code{r}" functions.
#'
#' @author Ben Goldstein, Perry de Valpine, and Daniel Turek
#'
#' @details
#'
#' These nimbleFunctions provide distributions that can be used directly in R or
#' in \code{nimble} hierarchical models (via \code{\link[nimble]{nimbleCode}}
#' and \code{\link[nimble]{nimbleModel}}).
#'
#' The letters following the 'dCJS_' indicate whether survival and/or
#' capture probabilities, in that order, are scalar (s, meaning the
#' probability applies to every \code{x[t]}) or vector (v, meaning the
#' probability is a vector aligned with \code{x}).  When
#' \code{probCapture} and/or \code{probSurvive} is a vector, they must
#' be the same length as \code{x}.
#'
#' It is important to use the time indexing correctly for survival.
#' \code{probSurvive[t]} is the survival probabilty from time
#' \code{t-1} to time \code{t}.  Time indexing for detection is more
#' obvious: \code{probDetect[t]} is the detection probability at time
#' \code{t}.
#'
#' When called from R, the \code{len} argument to \code{dCJS_**} is not
#' necessary. It will default to the length of \code{x}.  When used in
#' \code{nimble} model code (via \code{nimbleCode}), \code{len} must be provided
#' (even though it may seem redundant).
#'
#' For more explanation, see
#' \href{../doc/Introduction_to_nimbleEcology.html}{package vignette} (or
#' \code{vignette("Introduction_to_nimbleEcology")}).
#'
#' Compared to writing \code{nimble} models with a discrete latent state for
#' true alive/dead status at each time and a separate scalar datum for each observation, use
#' of these distributions allows one to directly sum (marginalize) over the
#' discrete latent states and calculate the probability of the detection history for one individual jointly.
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
#' \code{captures[i, 1:T] ~ dCSJ_ss(survive, capture, T)}
#'
#' declares a vector node, \code{captures[i, 1:T]}, (detection history for individual
#' \code{i},  for example) that follows a CJS distribution
#' with scalar survival probability \code{survive} and scalar capture probability \code{capture}
#' (assuming \code{survive} and \code{capture} are defined elsewhere in the model).
#'
#' This will invoke (something like) the following call to \code{dCJS_ss} when \code{nimble} uses the
#' model such as for MCMC:
#'
#' \code{dCJS_ss(captures[i, 1:T], survive, capture, len = T, log = TRUE)}
#'
#' If an algorithm using a \code{nimble} model with this declaration
#' needs to generate a random draw for \code{captures[i, 1:T]}, it
#' will make a similar invocation of \code{rCJS_ss}, with \code{n = 1}.
#'
#' If both survival and capture probabilities are time-dependent, use
#'
#' \code{captures[i,1:T] ~ dCSJ_vv(survive[1:T], capture[1:T], T)}
#'
#' and so on for each combination of time-dependent and time-independent parameters.
#'
#' @return
#'
#' For \code{dCJS_**}: the probability (or likelihood) or log probability of observation vector \code{x}.
#'
#' For \code{rCJS_**}: a simulated capture history, \code{x}.
#'
#' @references D. Turek, P. de Valpine and C. J. Paciorek. 2016. Efficient Markov chain Monte
#' Carlo sampling for hierarchical hidden Markov models. Environmental and Ecological Statistics
#' 23:549–564. DOI 10.1007/s10651-016-0353-z
#'
#' @seealso For multi-state or multi-event capture-recapture models, see \code{\link{dHMM}} or \code{\link{dDHMM}}.
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
#'   x[1:4] ~ dCJS_ss(probSurvive, probCapture, len = 4)
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
      if (!is.na(x[t])) {
        if (x[t] == 1) {
          ## ProbThisObs = P(x(t) | x(1)...x(t-1))
          probThisObs <- probAlive * probCapture
          probAliveGivenHistory <- 1
        } else {
          probAliveNotSeen <- probAlive * (1 - probCapture)
          probThisObs <- probAliveNotSeen + (1 - probAlive)
          probAliveGivenHistory <- probAliveNotSeen / probThisObs
        }
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
      if (!is.na(x[t])) {
        if (x[t] == 1) {
          ## ProbThisObs = P(x(t) | x(1)...x(t-1))
          probThisObs <- probAlive * probCapture[t]
          probAliveGivenHistory <- 1
        } else {
          probAliveNotSeen <- probAlive * (1 - probCapture[t])
          probThisObs <- probAliveNotSeen + (1 - probAlive)
          probAliveGivenHistory <- probAliveNotSeen / probThisObs
        }
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
      if (!is.na(x[t])) {
        if (x[t] == 1) {
          ## ProbThisObs = P(x(t) | x(1)...x(t-1))
          probThisObs <- probAlive * probCapture
          probAliveGivenHistory <- 1
        } else {
          probAliveNotSeen <- probAlive * (1 - probCapture)
          probThisObs <- probAliveNotSeen + (1 - probAlive)
          probAliveGivenHistory <- probAliveNotSeen / probThisObs
        }
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
      if (!is.na(x[t])) {
        if (x[t] == 1) {
          ## ProbThisObs = P(x(t) | x(1)...x(t-1))
          probThisObs <- probAlive * probCapture[t]
          probAliveGivenHistory <- 1
        } else {
          probAliveNotSeen <- probAlive * (1 - probCapture[t])
          probThisObs <- probAliveNotSeen + (1 - probAlive)
          probAliveGivenHistory <- probAliveNotSeen / probThisObs
        }
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
                 probSurvive = double(),
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
                 probCapture = double(),
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
    ans <- numeric(length = len, init = FALSE)
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
    types = c('value = double(1)', 'probSurvive = double()', 'probCapture = double()', 'len = double()'),
    pqAvail = FALSE))
  )

registerDistributions(list(
  dCJS_sv = list(
    BUGSdist = "dCJS_sv(probSurvive, probCapture, len)",
    Rdist = "dCJS_sv(probSurvive, probCapture, len = 0)",
    discrete = TRUE,
    types = c('value = double(1)', 'probSurvive = double()', 'probCapture = double(1)', 'len = double()'),
    pqAvail = FALSE))
  )

registerDistributions(list(
  dCJS_vs = list(
    BUGSdist = "dCJS_vs(probSurvive, probCapture, len)",
    Rdist = "dCJS_vs(probSurvive, probCapture, len = 0)",
    discrete = TRUE,
    types = c('value = double(1)', 'probSurvive = double(1)', 'probCapture = double()', 'len = double()'),
    pqAvail = FALSE))
  )

registerDistributions(list(
  dCJS_vv = list(
    BUGSdist = "dCJS_vv(probSurvive, probCapture, len)",
    Rdist = "dCJS_vv(probSurvive, probCapture, len = 0)",
    discrete = TRUE,
    types = c('value = double(1)', 'probSurvive = double(1)', 'probCapture = double(1)', 'len = double()'),
    pqAvail = FALSE))
)

