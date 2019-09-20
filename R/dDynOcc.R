# dDynOcc
#' Dynamic occupancy distribution for use in NIMBLE models
#'
#' \code{dDynOcc_**} provides dynamic occupancy model distributions for NIMBLE models.
#'
#' Dynamic occupancy models
#' model occurrence at a series of sites over many replicate timesteps. The likelihood of an observation in site s
#' at time t depends on the state of the site at time t-1, the transitition probability of persistence \code{probPersist[t]}
#' or colonization \code{probColonize[t]} and a detection probability code{p[s, t]}.
#'
#' The pair of letters following the 'dDynOcc_' indicates whether the probabilities of persistence
#' and colonization are a constant scalar (s) or time-indexed vector (v). For example, dOcc_sv takes scalar
#' persistence probability probPersist with a vector of colonization probabilities probColonize.
#'
#' Compared to writing NIMBLE models with a discrete latent state for true occupancy status and
#' a separate scalar datum for each observation,
#' use of these distributions allows
#' one to directly sum over the discrete latent state and calculate the probability of
#' all observations from one site jointly.
#'
#' @name dDynOcc
#' @aliases dDynOcc_ss dDynOcc_sv dDynOcc_vs dDynOcc_vv
#'
#' @param x detection/non-detection matrix of 0s (not detected) and 1s (detected). Each row contains repeat visits during one sampling period
#' @param init probability of occupancy in the first sampling period
#' @param probPersist persistence probability--probability an occupied cell remains occupied. 1-extinction probability. Scalar for \code{dDynOcc_s*}, vector for \code{dDynOcc_v*}
#' @param probColonize colonization probability. Probability that an unoccupied cell becomes occupied. \code{dDynOcc_*s}, vector for \code{dDynOcc*v}
#' @param p matrix of detection probabilities for each observation. Dimensions should match x
#' @param log TRUE (return log probability) or FALSE (return probability)
#' @param start a vector of the indices of the first observation in each time interval.
#' @param end a vector of the indices of the final observation in each time interval.
#' @param n a vector of length 2 indicating the dimensions of the data to be randomly generated
#'
#' @author Ben Goldstein and Perry de Valpine
#'
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
#' \code{detections[1:S, 1:T] ~ dDynOcc_ss(nrep, init = init_prob, probPersist = persistence_prob,
#' probColonize = colonization_prob, p = p[1:S, 1:T])}
#'
#' declares that the \code{detections[1:T]} vector follows a dynamic occupancy model distribution
#' with parameters as indicated, assuming all the parameters have been declared elsewhere in the model.
#'
#' If the colonization probabilities are time-dependent, one would use:
#'
#' \code{detections[1:T] ~ dDynOcc_sv(nrep, init = init_prob, probPersist = persistence_prob,
#' probColonize = colonization_prob[1:T], p = p[1:S, 1:T])}
#'
#' @seealso For regular occupancy models, see documentation for dOcc.
#' @examples
#' \dontrun{
#' # Set up constants and initial values for defining the model
#'   x <- matrix(c(0,0,NA,0,
#'                 1,1,1,0,
#'                 0,0,0,0,
#'                 0,0,1,0,
#'                 0,0,0,NA), nrow = 4)
#'   start <- c(1,1,2,1)
#'   end <- c(5,5,5,4)
#'   init <- 0.7
#'   probPersist <- 0.5
#'   probColonize <- 0.2
#'   p <- 0.8
#'
#'
#' # Define code for a nimbleModel
#'  nc <- nimbleCode({
#'
#'    x[1:2, 1:5] ~ dDynOcc_vv(nrep[1:2], init, probPersist[1:2], probColonize[1:2], p[1:2,1:5])
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
#' DynOcc_model <- nimbleModel(nc, data = list(x = dat, nrep = nrep),
#'                             inits = list(p = p, probPersist = probPersist,
#'                                          init = init, probColonize = probColonize))
#'
#' # Calculate log probability of data from the model
#' DynOcc_model$calculate()
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
    if (length(probPersist) != dim(x)[1] - 1) stop("Length of probPersist vector does not match length of data.")
    if (length(probColonize) != dim(x)[1] - 1) stop("Length of probColonize vector does not match length of data.")
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

    if (length(probPersist) != dim(x)[1] - 1) stop("Length of probPersist vector does not match length of data.")
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
    if (length(probColonize) != dim(x)[1] - 1) stop("Length of probColonize vector does not match length of data.")
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


registerDistributions(list(
    dDynOcc_vvm = list(
        BUGSdist = "dDynOcc_vvm(init, probPersist, probColonize, p, start, end)",
        Rdist = "dDynOcc_vvm(init, probPersist, probColonize, p, start, end)",
        types = c('value = double(2)',
                  'init = double()',
                  'probPersist = double(1)',
                  'probColonize = double(1)',
                  'p = double(2)',
                  'start = double(1)',
                  'end = double(1)'),
        mixedSizes = TRUE)))
registerDistributions(list(
    dDynOcc_vsm = list(
        BUGSdist = "dDynOcc_vsm(init, probPersist, probColonize, p, start, end)",
        Rdist = "dDynOcc_vsm(init, probPersist, probColonize, p, start, end)",
        types = c('value = double(2)',
                  'init = double()',
                  'probPersist = double(1)',
                  'probColonize = double()',
                  'p = double(2)',
                  'start = double(1)',
                  'end = double(1)'),
        mixedSizes = TRUE)))

registerDistributions(list(
    dDynOcc_svm = list(
        BUGSdist = "dDynOcc_svm(init, probPersist, probColonize, p, start, end)",
        Rdist = "dDynOcc_svm(init, probPersist, probColonize, p, start, end)",
        types = c('value = double(2)',
                  'init = double()',
                  'probPersist = double()',
                  'probColonize = double(1)',
                  'p = double(2)',
                  'start = double(1)',
                  'end = double(1)'),
        mixedSizes = TRUE)))

registerDistributions(list(
    dDynOcc_ssm = list(
        BUGSdist = "dDynOcc_ssm(init, probPersist, probColonize, p, start, end)",
        Rdist = "dDynOcc_ssm(init, probPersist, probColonize, p, start, end)",
        types = c('value = double(2)',
                  'init = double()',
                  'probPersist = double()',
                  'probColonize = double()',
                  'p = double(2)',
                  'start = double(1)',
                  'end = double(1)'),
        mixedSizes = TRUE)))




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
    if (length(probPersist) != dim(x)[1] - 1) stop("Length of probPersist vector does not match length of data.")
    if (length(probColonize) != dim(x)[1] - 1) stop("Length of probColonize vector does not match length of data.")
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

    if (length(probPersist) != dim(x)[1] - 1) stop("Length of probPersist vector does not match length of data.")
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
    if (length(probColonize) != dim(x)[1] - 1) stop("Length of probColonize vector does not match length of data.")
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


registerDistributions(list(
    dDynOcc_vvv = list(
        BUGSdist = "dDynOcc_vvv(init, probPersist, probColonize, p, start, end)",
        Rdist = "dDynOcc_vvv(init, probPersist, probColonize, p, start, end)",
        types = c('value = double(2)',
                  'init = double()',
                  'probPersist = double(1)',
                  'probColonize = double(1)',
                  'p = double(1)',
                  'start = double(1)',
                  'end = double(1)'),
        mixedSizes = TRUE)))
registerDistributions(list(
    dDynOcc_vsv = list(
        BUGSdist = "dDynOcc_vsv(init, probPersist, probColonize, p, start, end)",
        Rdist = "dDynOcc_vsv(init, probPersist, probColonize, p, start, end)",
        types = c('value = double(2)',
                  'init = double()',
                  'probPersist = double(1)',
                  'probColonize = double()',
                  'p = double(1)',
                  'start = double(1)',
                  'end = double(1)'),
        mixedSizes = TRUE)))

registerDistributions(list(
    dDynOcc_svv = list(
        BUGSdist = "dDynOcc_svv(init, probPersist, probColonize, p, start, end)",
        Rdist = "dDynOcc_svv(init, probPersist, probColonize, p, start, end)",
        types = c('value = double(2)',
                  'init = double()',
                  'probPersist = double()',
                  'probColonize = double(1)',
                  'p = double(1)',
                  'start = double(1)',
                  'end = double(1)'),
        mixedSizes = TRUE)))

registerDistributions(list(
    dDynOcc_ssv = list(
        BUGSdist = "dDynOcc_ssv(init, probPersist, probColonize, p, start, end)",
        Rdist = "dDynOcc_ssv(init, probPersist, probColonize, p, start, end)",
        types = c('value = double(2)',
                  'init = double()',
                  'probPersist = double()',
                  'probColonize = double()',
                  'p = double(1)',
                  'start = double(1)',
                  'end = double(1)'),
        mixedSizes = TRUE)))




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
    if (length(probPersist) != dim(x)[1] - 1) stop("Length of probPersist vector does not match length of data.")
    if (length(probColonize) != dim(x)[1] - 1) stop("Length of probColonize vector does not match length of data.")

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

    if (length(probPersist) != dim(x)[1] - 1) stop("Length of probPersist vector does not match length of data.")

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
    if (length(probColonize) != dim(x)[1] - 1) stop("Length of probColonize vector does not match length of data.")

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


registerDistributions(list(
    dDynOcc_vvs = list(
        BUGSdist = "dDynOcc_vvs(init, probPersist, probColonize, p, start, end)",
        Rdist = "dDynOcc_vvs(init, probPersist, probColonize, p, start, end)",
        types = c('value = double(2)',
                  'init = double()',
                  'probPersist = double(1)',
                  'probColonize = double(1)',
                  'p = double()',
                  'start = double(1)',
                  'end = double(1)'),
        mixedSizes = TRUE)))
registerDistributions(list(
    dDynOcc_vss = list(
        BUGSdist = "dDynOcc_vss(init, probPersist, probColonize, p, start, end)",
        Rdist = "dDynOcc_vss(init, probPersist, probColonize, p, start, end)",
        types = c('value = double(2)',
                  'init = double()',
                  'probPersist = double(1)',
                  'probColonize = double()',
                  'p = double()',
                  'start = double(1)',
                  'end = double(1)'),
        mixedSizes = TRUE)))

registerDistributions(list(
    dDynOcc_svs = list(
        BUGSdist = "dDynOcc_svs(init, probPersist, probColonize, p, start, end)",
        Rdist = "dDynOcc_svs(init, probPersist, probColonize, p, start, end)",
        types = c('value = double(2)',
                  'init = double()',
                  'probPersist = double()',
                  'probColonize = double(1)',
                  'p = double()',
                  'start = double(1)',
                  'end = double(1)'),
        mixedSizes = TRUE)))

registerDistributions(list(
    dDynOcc_sss = list(
        BUGSdist = "dDynOcc_sss(init, probPersist, probColonize, p, start, end)",
        Rdist = "dDynOcc_sss(init, probPersist, probColonize, p, start, end)",
        types = c('value = double(2)',
                  'init = double()',
                  'probPersist = double()',
                  'probColonize = double()',
                  'p = double()',
                  'start = double(1)',
                  'end = double(1)'),
        mixedSizes = TRUE)))

