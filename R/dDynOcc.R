# dDynOcc
#' Dynamic occupancy distribution for use in \code{nimble} models
#' \code{dDynOcc_**} and \code{rDynOcc_**} provide dynamic occupancy
#' model distributions that can be used directly from R or in \code{nimble}
#' models.
#'
#' @name dDynOcc
#' @aliases dDynOcc_sss dDynOcc_svs dDynOcc_vss dDynOcc_vvs
#' dDynOcc_ssv dDynOcc_svv dDynOcc_vsv dDynOcc_vvv
#' dDynOcc_ssm dDynOcc_svm dDynOcc_vsm dDynOcc_vvm
#' @author Ben Goldstein, Perry de Valpine and Lauren Ponisio
#'
#' @param x detection/non-detection matrix of 0s (not detected) and 1s
#'     (detected). Rows represent primary sampling occasions (e.g. different
#'     seasons). Columns are secondary sampling locations (e.g. replicate
#'     visits within a season) that may be different for each row
#' @param init probability of occupancy in the first sampling period
#' @param probPersist persistence probability--probability an occupied
#'     cell remains occupied. 1-extinction probability. Scalar for
#'     \code{dDynOcc_s**}, vector for \code{dDynOcc_v**}. If vector,
#'     should have length dim(x)[1] - 1 since no transition occurs
#'     after the last observation
#' @param probColonize colonization probability. Probability that
#'     an unoccupied cell becomes occupied. Scalar for \code{dDynOcc_*s*},
#'     vector for \code{dDynOcc_*v*}. If vector, should have length
#'     dim(x)[1] - 1 since no transition occurs after the last observation
#' @param p Detection probabilities. Scalar for \code{dDynOcc_**s},
#'     vector for \code{dDynOcc_**v}, matrix for \code{dDynOcc_**m}.
#'     If a matrix, dimensions should match x
#' @param log TRUE (return log probability) or FALSE (return probability)
#' @param start indicates the column number of the first observation in each
#'     row of x. A vector of length dim(x)[1]. This allows for different time
#'     periods to have different numbers of sampling occasions
#' @param end indicates the column number of the last observation in each
#'     row of x. A vector of length dim(x)[1]. This allows for different time
#'     periods to have different numbers of sampling occasions
#' @param n number of random draws, each returning a matrix of dimension
#'     \code{c(min(start), max(end))}. Currently only \code{n = 1} is supported,
#'     but the argument exists for standardization of "\code{r}" functions
#'
#' @details
#'
#' These nimbleFunctions provide distributions that can be used directly in R or
#' in \code{nimble} hierarchical models (via \code{\link[nimble]{nimbleCode}}
#' and \code{\link[nimble]{nimbleModel}}).
#'
#' The probability (or likelihood) of observation \code{x[t, o]} depends on
#' the occupancy status of the site at time t-1, the transitition
#' probability of persistence (\code{probPersist} or \code{probPersist[t]}),
#' colonization (\code{probColonize} or \code{probColonize[t]}), and a
#' detection probability (\code{p}, \code{p[t]}, or \code{p[t, o]}).
#'
#' The first two letters following the 'dDynOcc_' indicate whether the
#' probabilities of persistence and colonization are a constant scalar (s)
#' or time-indexed vector (v). For example, \code{dDynOcc_svm} takes scalar
#' persistence probability \code{probPersist} with a vector of colonization
#' probabilities \code{probColonize[1:(T-1)]}.
#'
#' When vectors, \code{probColonize} and \code{probPersist} may be of any
#' length greater than \code{length(x) - 1}. Only the first \code{length(x) - 1}
#' indices are used, each corresponding to the transition from time t to t+1
#' (e.g. \code{probColonize[2]} describes the transition probability from
#' t = 2 to t = 3). All extra values are ignored. This is to make it easier to
#' use one distribution for many sites, some requiring probabilities of length 1.
#'
#' The third letter in the suffix indicates whether the detection probability
#' is a constant (scalar), time-dependent (vector), or both time-dependent and
#' dependent on observation occasion (matrix). For example, \code{dDynOcc_svm}
#' takes a matrix of detection probabilities \code{p[1:T, 1:O]}.
#'
#' The arguments \code{start} and \code{end} allow different time periods to
#' contain different numbers of sampling events. Suppose you have observations
#' for samples in three seasons; in the first two seasons, there are four
#' observations, but in the third, there are only three. The \code{start}
#' and \code{end} could be provided as \code{start = c(1,1,1)} and
#' \code{end = c(4,4,3)}. In this case, the value of \code{x[4,4]} would
#' be ignored.
#'
#' For more explanation, see
#' \href{../doc/Introduction_to_nimbleEcology.html}{package vignette} (or
#' \code{vignette("Introduction_to_nimbleEcology")}).
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
#' \code{detections[1:T, 1:O] ~ dDynOcc_ssm(init,
#' probPersist = persistence_prob,
#' probColonize = colonization_prob, p = p[1:T, 1:O],
#' start = start[1:T], end = end[1:T])}
#'
#' declares that the \code{detections[1:T]} vector follows a dynamic occupancy
#' model distribution with parameters as indicated, assuming all the parameters
#' have been declared elsewhere in the model. This
#' will invoke (something like) the following call to \code{dDynOcc_ssm} when
#' \code{nimble} uses the model such as for MCMC:
#'
#' \code{dDynOcc_ssm(detections[1:T, 1:O], init,
#' probPersist = persistence_prob,
#' probColonize = colonization_prob, p = p[1:T, 1:O],
#' start = start[1:T], end = end[1:T], log = TRUE)}
#'
#' If an algorithm using a \code{nimble} model with this declaration
#' needs to generate a random draw for \code{detections[1:T, 1:O]}, it
#' will make a similar invocation of \code{rDynOcc_svm}, with \code{n = 1}.
#'
#' If the colonization probabilities are time-dependent, one would use:
#'
#' \code{detections[1:T] ~ dDynOcc_svm(nrep, init = init_prob,
#' probPersist = persistence_prob,
#' probColonize = colonization_prob[1:(T-1)], p = p[1:S, 1:T])}
#'
#' @return
#' For \code{dDynOcc_***}: the probability (or likelihood) or log probability
#' of observation vector \code{x}.
#' For \code{rDynOcc_***}: a simulated detection history, \code{x}.
#'
#' @seealso For basic occupancy models, see documentation for
#'   \code{\link{dOcc}}.
#' @examples
#' \donttest{
#' # Set up constants and initial values for defining the model
#'   x <- matrix(c(0,0,0,0,
#'                 1,1,1,0,
#'                 0,0,0,0,
#'                 0,0,1,0,
#'                 0,0,0,0), nrow = 4)
#'   start <- c(1,1,2,1,1)
#'   end <- c(5,5,5,4,5)
#'   init <- 0.7
#'   probPersist <- 0.5
#'   probColonize <- 0.2
#'   p <- matrix(rep(0.5, 20), nrow = 4)
#'
#'
#' # Define code for a nimbleModel
#'  nc <- nimbleCode({
#'
#'    x[1:2, 1:5] ~ dDynOcc_vvm(init,
#'      probPersist[1:2], probColonize[1:2], p[1:2,1:5],
#'      start = start[1:4], end = end[1:4])
#'
#'    init ~ dunif(0,1)
#'
#'    for (i in 1:2) {
#'      probPersist[i] ~ dunif(0,1)
#'      probColonize[i] ~ dunif(0,1)
#'    }
#'
#'    for (i in 1:2) {
#'      for (j in 1:5) {
#'        p[i,j] ~ dunif(0,1)
#'      }
#'    }
#'  })
#'
#' # Build the model, providing data and initial values
#' DynOcc_model <- nimbleModel(nc, data = list(x = x),
#'                             constants = list(start = start, end = end),
#'                             inits = list(p = p, probPersist = probPersist,
#'                                          init = init, probColonize = probColonize))
#'
#' # Calculate log probability of data from the model
#' DynOcc_model$calculate("x")
#' # Use the model for a variety of other purposes...
#' }

NULL
#' @rdname dDynOcc
#' @export
dDynOcc_vvm <- nimbleFunction(
  ## DynamicOccupancy removes the z's and muZ's from the model and computes
  ## the probability of all reps over all years for one site, species
  run = function(x = double(2),
                 init = double(),
                 probPersist = double(1),
                 probColonize = double(1),
                 p = double(2),
                 start = double(1),
                 end = double(1),
                 log = double(0, default = 0)) {
    if (length(probPersist) < dim(x)[1] - 1) stop("Length of probPersist vector must be at least length(x) - 1.")
    if (length(probColonize) < dim(x)[1] - 1) stop("Length of probColonize vector must be at least length(x) - 1.")
    if (dim(p)[1] != dim(x)[1]) stop("Dimension mismatch between x and p matrices.")
    if (dim(p)[2] != dim(x)[2]) stop("Dimension mismatch between x and p matrices.")

    ## x is a year by rep matix
    ProbOccNextTime <- init
    ll <- 0
    nyears <- dim(x)[1]
    if (nyears >= 1) {
      for (t in 1:nyears) {
        if (end[t] - start[t] + 1 > 0) {
          numObs <- sum(x[t, start[t]:end[t]])
          if (is.na(numObs)) numObs <- 0
          if (numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,start[t]:end[t]],
                             size = 1, prob = p[t,start[t]:end[t]], log = 1)))
          ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
          ProbCount <- ProbOccAndCount + ProbUnoccAndCount
          ProbOccGivenCount <- ProbOccAndCount / ProbCount
          ll <- ll + log(ProbCount)
          if (t < nyears)
              ProbOccNextTime <- ProbOccGivenCount * probPersist[t] +
                  (1 - ProbOccGivenCount) * probColonize[t]
        } else {
          ## If there were no observations in year t,
          ## simply propagate probability of occupancy forward
          if (t < nyears)
              ProbOccNextTime <- ProbOccNextTime *  probPersist[t] +
                  (1 - ProbOccNextTime) * probColonize[t]
        }
      }
    }
    if (log) return(ll)
    else return(exp(ll))
    returnType(double(0))
  }
)

#' @rdname dDynOcc
#' @export
dDynOcc_vsm <- nimbleFunction(
  ## DynamicOccupancy removes the z's and muZ's from the model and computes
  ## the probability of all reps over all years for one site, species
  run = function(x = double(2),
                 init = double(),
                 probPersist = double(1),
                 probColonize = double(),
                 p = double(2),
                 start = double(1),
                 end = double(1),
                 log = double(0, default = 0)) {

    if (length(probPersist) < dim(x)[1] - 1) stop("Length of probPersist vector must be at least length(x) - 1.")
    if (dim(p)[1] != dim(x)[1]) stop("Dimension mismatch between x and p matrices.")
    if (dim(p)[2] != dim(x)[2]) stop("Dimension mismatch between x and p matrices.")

    ## x is a year by rep matix
    ProbOccNextTime <- init
    ll <- 0
    nyears <- dim(x)[1]
    if (nyears >= 1) {
      for (t in 1:nyears) {
        if (end[t] - start[t] + 1 > 0) {
          numObs <- sum(x[t,start[t]:end[t]])
          if (numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,start[t]:end[t]],
                             size = 1, prob = p[t,start[t]:end[t]], log = 1)))
          ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
          ProbCount <- ProbOccAndCount + ProbUnoccAndCount
          ProbOccGivenCount <- ProbOccAndCount / ProbCount
          ll <- ll + log(ProbCount)
          if (t < nyears)
              ProbOccNextTime <- ProbOccGivenCount * probPersist[t] +
                  (1 - ProbOccGivenCount) * probColonize
        } else {
          ## If there were no observations in year t,
          ## simply propagate probability of occupancy forward
          if (t < nyears)
              ProbOccNextTime <- ProbOccNextTime *  probPersist[t] +
                  (1 - ProbOccNextTime) * probColonize
        }
      }
    }
    if (log) return(ll)
    else return(exp(ll))
    returnType(double(0))
  }
)

#' @rdname dDynOcc
#' @export
dDynOcc_svm <- nimbleFunction(
  ## DynamicOccupancy removes the z's and muZ's from the model and computes
  ## the probability of all reps over all years for one site, species
  run = function(x = double(2),
                 init = double(),
                 probPersist = double(),
                 probColonize = double(1),
                 p = double(2),
                 start = double(1),
                 end = double(1),
                 log = double(0, default = 0)) {
    if (length(probColonize) < dim(x)[1] - 1) stop("Length of probColonize vector must be at least length(x) - 1.")
    if (dim(p)[1] != dim(x)[1]) stop("Dimension mismatch between x and p matrices.")
    if (dim(p)[2] != dim(x)[2]) stop("Dimension mismatch between x and p matrices.")

    ## x is a year by rep matix
    ProbOccNextTime <- init
    ll <- 0
    nyears <- dim(x)[1]
    if (nyears >= 1) {
      for (t in 1:nyears) {
        if (end[t] - start[t] + 1 > 0) {
          numObs <- sum(x[t,start[t]:end[t]])
          if (numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,start[t]:end[t]],
                             size = 1, prob = p[t,start[t]:end[t]], log = 1)))
          ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
          ProbCount <- ProbOccAndCount + ProbUnoccAndCount
          ProbOccGivenCount <- ProbOccAndCount / ProbCount
          ll <- ll + log(ProbCount)
          if (t < nyears)
              ProbOccNextTime <- ProbOccGivenCount * probPersist +
                  (1 - ProbOccGivenCount) * probColonize[t]
        } else {
          ## If there were no observations in year t,
          ## simply propagate probability of occupancy forward
          if (t < nyears)
              ProbOccNextTime <- ProbOccNextTime * probPersist +
                  (1 - ProbOccNextTime) * probColonize[t]
        }
      }
    }
    if (log) return(ll)
    else return(exp(ll))
    returnType(double(0))
  }
)

#' @rdname dDynOcc
#' @export
dDynOcc_ssm <- nimbleFunction(
  ## DynamicOccupancy removes the z's and muZ's from the model and computes
  ## the probability of all reps over all years for one site, species
  run = function(x = double(2),
                 init = double(),
                 probPersist = double(),
                 probColonize = double(),
                 p = double(2),
                 start = double(1),
                 end = double(1),
                 log = double(0, default = 0)) {
    # if (length(probColonize) != 1) stop("In dDynOcc_vs probColonize must be scalar")
    # if (length(probPersist) != 1) stop("In dDynOcc_vs probPersist must be scalar")
    if (dim(p)[1] != dim(x)[1]) stop("Dimension mismatch between x and p matrices.")
    if (dim(p)[2] != dim(x)[2]) stop("Dimension mismatch between x and p matrices.")

    ## x is a year by rep matix
    ProbOccNextTime <- init
    ll <- 0
    nyears <- dim(x)[1]
    if (nyears >= 1) {
      for (t in 1:nyears) {
        if (end[t] - start[t] + 1 > 0) {
          numObs <- sum(x[t,start[t]:end[t]])
          if (numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,start[t]:end[t]],
                             size = 1, prob = p[t,start[t]:end[t]], log = 1)))
          ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
          ProbCount <- ProbOccAndCount + ProbUnoccAndCount
          ProbOccGivenCount <- ProbOccAndCount / ProbCount
          ll <- ll + log(ProbCount)
          if (t < nyears)
              ProbOccNextTime <- ProbOccGivenCount * probPersist +
                  (1 - ProbOccGivenCount) * probColonize
        } else {
          ## If there were no observations in year t,
          ## simply propagate probability of occupancy forward
          if (t < nyears)
              ProbOccNextTime <- ProbOccNextTime * probPersist +
                  (1 - ProbOccNextTime) * probColonize
        }
      }
    }
    if (log) return(ll)
    else return(exp(ll))
    returnType(double(0))
  }
)

#' @rdname dDynOcc
#' @export
rDynOcc_vvm <- nimbleFunction(
  run = function(n = double(),
                 init = double(),
                 probPersist = double(1),
                 probColonize = double(1),
                 p = double(2),
                 start = double(1),
                 end = double(1)) {
    occupied <- rbinom(1, 1, init)
    val <- matrix(-1, nrow = dim(p)[1], ncol = dim(p)[2])
    val[1, start[1]:end[1]] <- occupied * rbinom(end[1] - start[1] + 1, 1, p[1,])

    for (t in 2:dim(p)[1]) {
      if (occupied == 1) {
        occupied <- rbinom(1, 1, probPersist[t - 1])
      } else {
        occupied <- rbinom(1, 1, probColonize[t - 1])
      }
      val[t, start[t]:end[t]] <- occupied * rbinom(end[t] - start[t] + 1, 1, p[t,])
    }

    return(val)
    returnType(double(2))
  }
)

#' @rdname dDynOcc
#' @export
rDynOcc_vsm <- nimbleFunction(
  run = function(n = double(0),
                 init = double(),
                 probPersist = double(1),
                 probColonize = double(),
                 p = double(2),
                 start = double(1),
                 end = double(1)) {
    occupied <- rbinom(1, 1, init)
    val <- matrix(-1, nrow = dim(p)[1], ncol = dim(p)[2])
    val[1, start[1]:end[1]] <- occupied * rbinom(end[1] - start[1] + 1, 1, p[1,])

    for (t in 2:dim(p)[1]) {
      if (occupied == 1) {
        occupied <- rbinom(1, 1, probPersist[t - 1])
      } else {
        occupied <- rbinom(1, 1, probColonize)
      }
      val[t, start[t]:end[t]] <- occupied * rbinom(end[t] - start[t] + 1, 1, p[t,])
    }

    return(val)
    returnType(double(2))
  }
)

#' @rdname dDynOcc
#' @export
rDynOcc_svm <- nimbleFunction(
  run = function(n = double(0),
                 init = double(),
                 probPersist = double(),
                 probColonize = double(1),
                 p = double(2),
                 start = double(1),
                 end = double(1)) {
    occupied <- rbinom(1, 1, init)
    val <- matrix(-1, nrow = dim(p)[1], ncol = dim(p)[2])
    val[1, start[1]:end[1]] <- occupied * rbinom(end[1] - start[1] + 1, 1, p[1,])

    for (t in 2:dim(p)[1]) {
      if (occupied == 1) {
        occupied <- rbinom(1, 1, probPersist)
      } else {
        occupied <- rbinom(1, 1, probColonize[t - 1])
      }
      val[t, start[t]:end[t]] <- occupied * rbinom(end[t] - start[t] + 1, 1, p[t,])
    }

    return(val)
    returnType(double(2))
  }
)

#' @rdname dDynOcc
#' @export
rDynOcc_ssm <- nimbleFunction(
  run = function(n = double(0),
                 init = double(),
                 probPersist = double(),
                 probColonize = double(),
                 p = double(2),
                 start = double(1),
                 end = double(1)) {
    occupied <- rbinom(1, 1, init)
    val <- matrix(-1, nrow = dim(p)[1], ncol = dim(p)[2])
    val[1, start[1]:end[1]] <- occupied * rbinom(end[1] - start[1] + 1, 1, p[1,])

    for (t in 2:dim(p)[1]) {
      if (occupied == 1) {
        occupied <- rbinom(1, 1, probPersist)
      } else {
        occupied <- rbinom(1, 1, probColonize)
      }
      val[t, start[t]:end[t]] <- occupied * rbinom(end[t] - start[t] + 1, 1, p[t,])
    }

    return(val)
    returnType(double(2))
  }
)

NULL
#' @rdname dDynOcc
#' @export
dDynOcc_vvv <- nimbleFunction(
  ## DynamicOccupancy removes the z's and muZ's from the model and computes
  ## the probability of all reps over all years for one site, species
  run = function(x = double(2),
                 init = double(),
                 probPersist = double(1),
                 probColonize = double(1),
                 p = double(1),
                 start = double(1),
                 end = double(1),
                 log = double(0, default = 0)) {
    if (length(probPersist) < dim(x)[1] - 1) stop("Length of probPersist vector must be at least length(x) - 1.")
    if (length(probColonize) < dim(x)[1] - 1) stop("Length of probColonize vector must be at least length(x) - 1.")
    if (dim(p)[1] != dim(x)[1]) stop("Dimension mismatch between x matrix and p vector.")

    ## x is a year by rep matix
    ProbOccNextTime <- init
    ll <- 0
    nyears <- dim(x)[1]
    if (nyears >= 1) {
      for (t in 1:nyears) {
        if (end[t] - start[t] + 1 > 0) {
          numObs <- sum(x[t,start[t]:end[t]])
          if (numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,start[t]:end[t]],
                             size = 1, prob = p[t], log = 1)))
          ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
          ProbCount <- ProbOccAndCount + ProbUnoccAndCount
          ProbOccGivenCount <- ProbOccAndCount / ProbCount
          ll <- ll + log(ProbCount)
          if (t < nyears)
              ProbOccNextTime <- ProbOccGivenCount * probPersist[t] +
                  (1 - ProbOccGivenCount) * probColonize[t]
        } else {
          ## If there were no observations in year t,
          ## simply propagate probability of occupancy forward
          if (t < nyears)
              ProbOccNextTime <- ProbOccNextTime *  probPersist[t] +
                  (1 - ProbOccNextTime) * probColonize[t]
        }
      }
    }
    if (log) return(ll)
    else return(exp(ll))
    returnType(double(0))
  }
)

#' @rdname dDynOcc
#' @export
dDynOcc_vsv <- nimbleFunction(
  ## DynamicOccupancy removes the z's and muZ's from the model and computes
  ## the probability of all reps over all years for one site, species
  run = function(x = double(2),
                 init = double(),
                 probPersist = double(1),
                 probColonize = double(),
                 p = double(1),
                 start = double(1),
                 end = double(1),
                 log = double(0, default = 0)) {

    if (length(probPersist) < dim(x)[1] - 1) stop("Length of probPersist vector must be at least length(x) - 1.")
    if (dim(p)[1] != dim(x)[1]) stop("Dimension mismatch between x matrix and p vector.")

    ## x is a year by rep matix
    ProbOccNextTime <- init
    ll <- 0
    nyears <- dim(x)[1]
    if (nyears >= 1) {
      for (t in 1:nyears) {
        if (end[t] - start[t] + 1 > 0) {
          numObs <- sum(x[t,start[t]:end[t]])
          if (numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,start[t]:end[t]],
                             size = 1, prob = p[t], log = 1)))
          ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
          ProbCount <- ProbOccAndCount + ProbUnoccAndCount
          ProbOccGivenCount <- ProbOccAndCount / ProbCount
          ll <- ll + log(ProbCount)
          if (t < nyears)
              ProbOccNextTime <- ProbOccGivenCount * probPersist[t] +
                  (1 - ProbOccGivenCount) * probColonize
        } else {
          ## If there were no observations in year t,
          ## simply propagate probability of occupancy forward
          if (t < nyears)
              ProbOccNextTime <- ProbOccNextTime *  probPersist[t] +
                  (1 - ProbOccNextTime) * probColonize
        }
      }
    }
    if (log) return(ll)
    else return(exp(ll))
    returnType(double(0))
  }
)

#' @rdname dDynOcc
#' @export
dDynOcc_svv <- nimbleFunction(
  ## DynamicOccupancy removes the z's and muZ's from the model and computes
  ## the probability of all reps over all years for one site, species
  run = function(x = double(2),
                 init = double(),
                 probPersist = double(),
                 probColonize = double(1),
                 p = double(1),
                 start = double(1),
                 end = double(1),
                 log = double(0, default = 0)) {
    if (length(probColonize) < dim(x)[1] - 1) stop("Length of probColonize vector must be at least length(x) - 1.")
    if (dim(p)[1] != dim(x)[1]) stop("Dimension mismatch between x matrix and p vector.")

    ## x is a year by rep matix
    ProbOccNextTime <- init
    ll <- 0
    nyears <- dim(x)[1]
    if (nyears >= 1) {
      for (t in 1:nyears) {
        if (end[t] - start[t] + 1 > 0) {
          numObs <- sum(x[t,start[t]:end[t]])
          if (numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,start[t]:end[t]],
                             size = 1, prob = p[t], log = 1)))
          ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
          ProbCount <- ProbOccAndCount + ProbUnoccAndCount
          ProbOccGivenCount <- ProbOccAndCount / ProbCount
          ll <- ll + log(ProbCount)
          if (t < nyears)
              ProbOccNextTime <- ProbOccGivenCount * probPersist +
                  (1 - ProbOccGivenCount) * probColonize[t]
        } else {
          ## If there were no observations in year t,
          ## simply propagate probability of occupancy forward
          if (t < nyears)
              ProbOccNextTime <- ProbOccNextTime * probPersist +
                  (1 - ProbOccNextTime) * probColonize[t]
        }
      }
    }
    if (log) return(ll)
    else return(exp(ll))
    returnType(double(0))
  }
)

#' @rdname dDynOcc
#' @export
dDynOcc_ssv <- nimbleFunction(
  ## DynamicOccupancy removes the z's and muZ's from the model and computes
  ## the probability of all reps over all years for one site, species
  run = function(x = double(2),
                 init = double(),
                 probPersist = double(),
                 probColonize = double(),
                 p = double(1),
                 start = double(1),
                 end = double(1),
                 log = double(0, default = 0)) {
    # if (length(probColonize) != 1) stop("In dDynOcc_vs probColonize must be scalar")
    # if (length(probPersist) != 1) stop("In dDynOcc_vs probPersist must be scalar")
    if (dim(p)[1] != dim(x)[1]) stop("Dimension mismatch between x matrix and p vector.")

    ## x is a year by rep matix
    ProbOccNextTime <- init
    ll <- 0
    nyears <- dim(x)[1]
    if (nyears >= 1) {
      for (t in 1:nyears) {
        if (end[t] - start[t] + 1 > 0) {
          numObs <- sum(x[t,start[t]:end[t]])
          if (numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,start[t]:end[t]],
                             size = 1, prob = p[t], log = 1)))
          ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
          ProbCount <- ProbOccAndCount + ProbUnoccAndCount
          ProbOccGivenCount <- ProbOccAndCount / ProbCount
          ll <- ll + log(ProbCount)
          if (t < nyears)
              ProbOccNextTime <- ProbOccGivenCount * probPersist +
                  (1 - ProbOccGivenCount) * probColonize
        } else {
          ## If there were no observations in year t,
          ## simply propagate probability of occupancy forward
          if (t < nyears)
              ProbOccNextTime <- ProbOccNextTime * probPersist +
                  (1 - ProbOccNextTime) * probColonize
        }
      }
    }
    if (log) return(ll)
    else return(exp(ll))
    returnType(double(0))
  }
)


#' @rdname dDynOcc
#' @export
rDynOcc_vvv <- nimbleFunction(
  run = function(n = double(),
                 init = double(),
                 probPersist = double(1),
                 probColonize = double(1),
                 p = double(1),
                 start = double(1),
                 end = double(1)) {
    occupied <- rbinom(1, 1, init)
    val <- matrix(-1, nrow = dim(p)[1], ncol = max(end))
    val[1, start[1]:end[1]] <- occupied * rbinom(end[1] - start[1] + 1, 1, p[1])

    for (t in 2:dim(p)[1]) {
      if (occupied == 1) {
        occupied <- rbinom(1, 1, probPersist[t - 1])
      } else {
        occupied <- rbinom(1, 1, probColonize[t - 1])
      }
      val[t, start[t]:end[t]] <- occupied * rbinom(end[t] - start[t] + 1, 1, p[t])
    }

    return(val)
    returnType(double(2))
  }
)

#' @rdname dDynOcc
#' @export
rDynOcc_vsv <- nimbleFunction(
  run = function(n = double(),
                 init = double(),
                 probPersist = double(1),
                 probColonize = double(),
                 p = double(1),
                 start = double(1),
                 end = double(1)) {
    occupied <- rbinom(1, 1, init)
    val <- matrix(-1, nrow = dim(p)[1], ncol = max(end))
    val[1, start[1]:end[1]] <- occupied * rbinom(end[1] - start[1] + 1, 1, p[1])

    for (t in 2:dim(p)[1]) {
      if (occupied == 1) {
        occupied <- rbinom(1, 1, probPersist[t - 1])
      } else {
        occupied <- rbinom(1, 1, probColonize)
      }
      val[t, start[t]:end[t]] <- occupied * rbinom(end[t] - start[t] + 1, 1, p[t])
    }

    return(val)
    returnType(double(2))
  }
)

#' @rdname dDynOcc
#' @export
rDynOcc_svv <- nimbleFunction(
  run = function(n = double(),
                 init = double(),
                 probPersist = double(),
                 probColonize = double(1),
                 p = double(1),
                 start = double(1),
                 end = double(1)) {
    occupied <- rbinom(1, 1, init)
    val <- matrix(-1, nrow = dim(p)[1], ncol = max(end))
    val[1, start[1]:end[1]] <- occupied * rbinom(end[1] - start[1] + 1, 1, p[1])

    for (t in 2:dim(p)[1]) {
      if (occupied == 1) {
        occupied <- rbinom(1, 1, probPersist)
      } else {
        occupied <- rbinom(1, 1, probColonize[t - 1])
      }
      val[t, start[t]:end[t]] <- occupied * rbinom(end[t] - start[t] + 1, 1, p[t])
    }

    return(val)
    returnType(double(2))
  }
)

#' @rdname dDynOcc
#' @export
rDynOcc_ssv <- nimbleFunction(
  run = function(n = double(),
                 init = double(),
                 probPersist = double(),
                 probColonize = double(),
                 p = double(1),
                 start = double(1),
                 end = double(1)) {
    occupied <- rbinom(1, 1, init)
    val <- matrix(-1, nrow = dim(p)[1], ncol = max(end))
    val[1, start[1]:end[1]] <- occupied * rbinom(end[1] - start[1] + 1, 1, p[1])

    for (t in 2:dim(p)[1]) {
      if (occupied == 1) {
        occupied <- rbinom(1, 1, probPersist)
      } else {
        occupied <- rbinom(1, 1, probColonize)
      }
      val[t, start[t]:end[t]] <- occupied * rbinom(end[t] - start[t] + 1, 1, p[t])
    }

    return(val)
    returnType(double(2))
  }
)


NULL
#' @rdname dDynOcc
#' @export
dDynOcc_vvs <- nimbleFunction(
  ## DynamicOccupancy removes the z's and muZ's from the model and computes
  ## the probability of all reps over all years for one site, species
  run = function(x = double(2),
                 init = double(),
                 probPersist = double(1),
                 probColonize = double(1),
                 p = double(),
                 start = double(1),
                 end = double(1),
                 log = double(0, default = 0)) {
    if (length(probPersist) < dim(x)[1] - 1) stop("Length of probPersist vector must be at least length(x) - 1.")
    if (length(probColonize) < dim(x)[1] - 1) stop("Length of probColonize vector must be at least length(x) - 1.")

    ## x is a year by rep matix
    ProbOccNextTime <- init
    ll <- 0
    nyears <- dim(x)[1]
    if (nyears >= 1) {
      for (t in 1:nyears) {
        if (end[t] - start[t] + 1 > 0) {
          numObs <- sum(x[t,start[t]:end[t]])
          if (numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,start[t]:end[t]],
                             size = 1, prob = p, log = 1)))
          ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
          ProbCount <- ProbOccAndCount + ProbUnoccAndCount
          ProbOccGivenCount <- ProbOccAndCount / ProbCount
          ll <- ll + log(ProbCount)
          if (t < nyears)
              ProbOccNextTime <- ProbOccGivenCount * probPersist[t] +
                  (1 - ProbOccGivenCount) * probColonize[t]
        } else {
          ## If there were no observations in year t,
          ## simply propagate probability of occupancy forward
          if (t < nyears)
              ProbOccNextTime <- ProbOccNextTime *  probPersist[t] +
                  (1 - ProbOccNextTime) * probColonize[t]
        }
      }
    }
    if (log) return(ll)
    else return(exp(ll))
    returnType(double(0))
  }
)

#' @rdname dDynOcc
#' @export
dDynOcc_vss <- nimbleFunction(
  ## DynamicOccupancy removes the z's and muZ's from the model and computes
  ## the probability of all reps over all years for one site, species
  run = function(x = double(2),
                 init = double(),
                 probPersist = double(1),
                 probColonize = double(),
                 p = double(),
                 start = double(1),
                 end = double(1),
                 log = double(0, default = 0)) {

    if (length(probPersist) < dim(x)[1] - 1) stop("Length of probPersist vector must be at least length(x) - 1.")

    ## x is a year by rep matix
    ProbOccNextTime <- init
    ll <- 0
    nyears <- dim(x)[1]
    if (nyears >= 1) {
      for (t in 1:nyears) {
        if (end[t] - start[t] + 1 > 0) {
          numObs <- sum(x[t,start[t]:end[t]])
          if (numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,start[t]:end[t]],
                             size = 1, prob = p, log = 1)))
          ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
          ProbCount <- ProbOccAndCount + ProbUnoccAndCount
          ProbOccGivenCount <- ProbOccAndCount / ProbCount
          ll <- ll + log(ProbCount)
          if (t < nyears)
              ProbOccNextTime <- ProbOccGivenCount * probPersist[t] +
                  (1 - ProbOccGivenCount) * probColonize
        } else {
          ## If there were no observations in year t,
          ## simply propagate probability of occupancy forward
          if (t < nyears)
              ProbOccNextTime <- ProbOccNextTime *  probPersist[t] +
                  (1 - ProbOccNextTime) * probColonize
        }
      }
    }
    if (log) return(ll)
    else return(exp(ll))
    returnType(double(0))
  }
)

#' @rdname dDynOcc
#' @export
dDynOcc_svs <- nimbleFunction(
  ## DynamicOccupancy removes the z's and muZ's from the model and computes
  ## the probability of all reps over all years for one site, species
  run = function(x = double(2),
                 init = double(),
                 probPersist = double(),
                 probColonize = double(1),
                 p = double(),
                 start = double(1),
                 end = double(1),
                 log = double(0, default = 0)) {
    if (length(probColonize) < dim(x)[1] - 1) stop("Length of probColonize vector must be at least length(x) - 1.")

    ## x is a year by rep matix
    ProbOccNextTime <- init
    ll <- 0
    nyears <- dim(x)[1]
    if (nyears >= 1) {
      for (t in 1:nyears) {
        if (end[t] - start[t] + 1 > 0) {
          numObs <- sum(x[t,start[t]:end[t]])
          if (numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,start[t]:end[t]],
                             size = 1, prob = p, log = 1)))
          ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
          ProbCount <- ProbOccAndCount + ProbUnoccAndCount
          ProbOccGivenCount <- ProbOccAndCount / ProbCount
          ll <- ll + log(ProbCount)
          if (t < nyears)
              ProbOccNextTime <- ProbOccGivenCount * probPersist +
                  (1 - ProbOccGivenCount) * probColonize[t]
        } else {
          ## If there were no observations in year t,
          ## simply propagate probability of occupancy forward
          if (t < nyears)
              ProbOccNextTime <- ProbOccNextTime * probPersist +
                  (1 - ProbOccNextTime) * probColonize[t]
        }
      }
    }
    if (log) return(ll)
    else return(exp(ll))
    returnType(double(0))
  }
)

#' @rdname dDynOcc
#' @export
dDynOcc_sss <- nimbleFunction(
  ## DynamicOccupancy removes the z's and muZ's from the model and computes
  ## the probability of all reps over all years for one site, species
  run = function(x = double(2),
                 init = double(),
                 probPersist = double(),
                 probColonize = double(),
                 p = double(),
                 start = double(1),
                 end = double(1),
                 log = double(0, default = 0)) {

    ## x is a year by rep matix
    ProbOccNextTime <- init
    ll <- 0
    nyears <- dim(x)[1]
    if (nyears >= 1) {
      for (t in 1:nyears) {
        if (end[t] - start[t] + 1 > 0) {
          numObs <- sum(x[t,start[t]:end[t]])
          if (numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,start[t]:end[t]],
                             size = 1, prob = p, log = 1)))
          ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
          ProbCount <- ProbOccAndCount + ProbUnoccAndCount
          ProbOccGivenCount <- ProbOccAndCount / ProbCount
          ll <- ll + log(ProbCount)
          if (t < nyears)
              ProbOccNextTime <- ProbOccGivenCount * probPersist +
                  (1 - ProbOccGivenCount) * probColonize
        } else {
          ## If there were no observations in year t,
          ## simply propagate probability of occupancy forward
          if (t < nyears)
              ProbOccNextTime <- ProbOccNextTime * probPersist +
                  (1 - ProbOccNextTime) * probColonize
        }
      }
    }
    if (log) return(ll)
    else return(exp(ll))
    returnType(double(0))
  }
)


#' @rdname dDynOcc
#' @export
rDynOcc_vvs <- nimbleFunction(
  run = function(n = double(),
                 init = double(),
                 probPersist = double(1),
                 probColonize = double(1),
                 p = double(),
                 start = double(1),
                 end = double(1)) {
    occupied <- rbinom(1, 1, init)
    val <- matrix(-1, nrow = length(end), ncol = max(end))
    val[1, start[1]:end[1]] <- occupied * rbinom(end[1] - start[1] + 1, 1, p)

    for (t in 2:length(end)) {
      if (occupied == 1) {
        occupied <- rbinom(1, 1, probPersist[t - 1])
      } else {
        occupied <- rbinom(1, 1, probColonize[t - 1])
      }
      val[t, start[t]:end[t]] <- occupied * rbinom(end[t] - start[t] + 1, 1, p)
    }

    return(val)
    returnType(double(2))
  }
)

#' @rdname dDynOcc
#' @export
rDynOcc_vss <- nimbleFunction(
  run = function(n = double(),
                 init = double(),
                 probPersist = double(1),
                 probColonize = double(),
                 p = double(),
                 start = double(1),
                 end = double(1)) {
    occupied <- rbinom(1, 1, init)
    val <- matrix(-1, nrow = length(end), ncol = max(end))
    val[1, start[1]:end[1]] <- occupied * rbinom(end[1] - start[1] + 1, 1, p)

    for (t in 2:length(end)) {
      if (occupied == 1) {
        occupied <- rbinom(1, 1, probPersist[t - 1])
      } else {
        occupied <- rbinom(1, 1, probColonize)
      }
      val[t, start[t]:end[t]] <- occupied * rbinom(end[t] - start[t] + 1, 1, p)
    }

    return(val)
    returnType(double(2))
  }
)

#' @rdname dDynOcc
#' @export
rDynOcc_svs <- nimbleFunction(
  run = function(n = double(),
                 init = double(),
                 probPersist = double(),
                 probColonize = double(1),
                 p = double(),
                 start = double(1),
                 end = double(1)) {
    occupied <- rbinom(1, 1, init)
    val <- matrix(-1, nrow = length(end), ncol = max(end))
    val[1, start[1]:end[1]] <- occupied * rbinom(end[1] - start[1] + 1, 1, p)

    for (t in 2:length(end)) {
      if (occupied == 1) {
        occupied <- rbinom(1, 1, probPersist)
      } else {
        occupied <- rbinom(1, 1, probColonize[t - 1])
      }
      val[t, start[t]:end[t]] <- occupied * rbinom(end[t] - start[t] + 1, 1, p)
    }

    return(val)
    returnType(double(2))
  }
)

#' @rdname dDynOcc
#' @export
rDynOcc_sss <- nimbleFunction(
  run = function(n = double(),
                 init = double(),
                 probPersist = double(),
                 probColonize = double(),
                 p = double(),
                 start = double(1),
                 end = double(1)) {
    occupied <- rbinom(1, 1, init)
    val <- matrix(-1, nrow = length(end), ncol = max(end))
    val[1, start[1]:end[1]] <- occupied * rbinom(end[1] - start[1] + 1, 1, p)

    for (t in 2:length(end)) {
      if (occupied == 1) {
        occupied <- rbinom(1, 1, probPersist)
      } else {
        occupied <- rbinom(1, 1, probColonize)
      }
      val[t, start[t]:end[t]] <- occupied * rbinom(end[t] - start[t] + 1, 1, p)
    }

    return(val)
    returnType(double(2))
  }
)


