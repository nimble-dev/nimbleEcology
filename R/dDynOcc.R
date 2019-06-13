# dDynOcc
# To be filled in with dynamic occupancy model distribution(s).
#' Dynamic occupancy distribution for use in NIMBLE models
#'
#' @aliases dDynOcc_ss dDynOcc_sv dDynOcc_vs dDynOcc_vv
#'
#' \code{dDynOcc_**} provides dynamic occupancy model distributions for NIMBLE models.
#' Dynamic occupancy models
#' The pair of letters following the 'dOcc_' indicates whether the probabilities of persistence
#' and colonization are scalar (s, uniform for all ) or vector (v). For example, dOcc_sc takes scalar
#' occupancy probability with a vector of detection probabilities.
#'
#' Compared to writing NIMBLE models with a discrete latent state for true occupancy status and
#' a separate scalar datum for each observation,
#' use of these distributions allows
#' one to directly sum over the discrete latent state and calculate the probability of
#' all observations from one site jointly.
#' @name dDynOcc
#'
#' @export
#'
#' @param x detection/non-detection matrix of 0s (not detected) and 1s (detected). Each row contains repeat visits during one sampling period
#' @param nrep a vector of the number of observations per sampling occasion
#' @param psi1 probability of occupancy in the first sampling period
#' @param phi persistence probability--probability an occupied cell remains occupied. 1-extinction probability. Scalar for \code{dOcc_s\*}, vector for \code{dOcc_v\*}
#' @param gamma colonization probability. Probability that an unoccupied cell becomes occupied. \code{dOcc_\*s}, vector for \code{dOcc_\*v}
#' @param p matrix of detection probabilities for each observation. Dimensions should match x
#' @param log TRUE (return log probability) or FALSE (return probability)
#'
#' @author Ben Goldstein and Perry de Valpine
#'
#' @export
dDynOcc_vv <- nimbleFunction(
  ## DynamicOccupancy removes the z's and muZ's from the model and computes
  ## the probability of all reps over all years for one site, species
  run = function(x = double(2),
                 nrep = double(1),
                 psi1 = double(0),
                 phi = double(1),
                 gamma = double(1),
                 p = double(2),
                 log = double(0, default = 0)) {
    if (length(phi) != dim(x)[1]) stop("Length of phi vector does not match length of data.")
    if (length(gamma) != dim(x)[1]) stop("Length of gamma vector does not match length of data.")
    if (dim(p)[1] != dim(x)[1]) stop("Dimension mismatch between x and p matrices.")
    if (dim(p)[2] != dim(x)[2]) stop("Dimension mismatch between x and p matrices.")

    ## x is a year by rep matix
    ProbOccNextTime <- psi1
    ll <- 0
    nyears <- dim(x)[1]
    if(nyears >= 1) {
      for(t in 1:nyears) {
        if(nrep[t] > 0) {
          numObs <- sum(x[t,1:nrep[t]])
          if(numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,1:nrep[t]],
                             size = 1, p = p[t,1:nrep[t]], log = 1)))
          ProbUnoccAndCount <- (1-ProbOccNextTime) * (numObs == 0)
          ProbCount <- ProbOccAndCount + ProbUnoccAndCount
          ProbOccGivenCount <- ProbOccAndCount / ProbCount
          ll <- ll + log(ProbCount)
          if(t < nyears)
              ProbOccNextTime <- ProbOccGivenCount * phi[t] +
                  (1-ProbOccGivenCount) * gamma[t]
        } else {
          ## If there were no observations in year t,
          ## simply propagate probability of occupancy forward
          if (t < nyears)
              ProbOccNextTime <- ProbOccNextTime *  phi[t] +
                  (1-ProbOccNextTime) * gamma[t]
        }
      }
    }
    if(log) return(ll)
    else return(exp(ll))
    returnType(double(0))
  }
)

dDynOcc_vs <- nimbleFunction(
  ## DynamicOccupancy removes the z's and muZ's from the model and computes
  ## the probability of all reps over all years for one site, species
  run = function(x = double(2),
                 nrep = double(1),
                 psi1 = double(0),
                 phi = double(1),
                 gamma = double(0),
                 p = double(2),
                 log = double(0, default = 0)) {

    if (length(phi) != dim(x)[1]) stop("Length of phi vector does not match length of data.")
    if (dim(p)[1] != dim(x)[1]) stop("Dimension mismatch between x and p matrices.")
    if (dim(p)[2] != dim(x)[2]) stop("Dimension mismatch between x and p matrices.")

    ## x is a year by rep matix
    ProbOccNextTime <- psi1
    ll <- 0
    nyears <- dim(x)[1]
    if(nyears >= 1) {
      for(t in 1:nyears) {
        if(nrep[t] > 0) {
          numObs <- sum(x[t,1:nrep[t]])
          if(numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,1:nrep[t]],
                             size = 1, p = p[t,1:nrep[t]], log = 1)))
          ProbUnoccAndCount <- (1-ProbOccNextTime) * (numObs == 0)
          ProbCount <- ProbOccAndCount + ProbUnoccAndCount
          ProbOccGivenCount <- ProbOccAndCount / ProbCount
          ll <- ll + log(ProbCount)
          if(t < nyears)
              ProbOccNextTime <- ProbOccGivenCount * phi[t] +
                  (1-ProbOccGivenCount) * gamma
        } else {
          ## If there were no observations in year t,
          ## simply propagate probability of occupancy forward
          if (t < nyears)
              ProbOccNextTime <- ProbOccNextTime *  phi[t] +
                  (1-ProbOccNextTime) * gamma
        }
      }
    }
    if(log) return(ll)
    else return(exp(ll))
    returnType(double(0))
  }
)


dDynOcc_sv <- nimbleFunction(
  ## DynamicOccupancy removes the z's and muZ's from the model and computes
  ## the probability of all reps over all years for one site, species
  run = function(x = double(2),
                 nrep = double(1),
                 psi1 = double(0),
                 phi = double(0),
                 gamma = double(1),
                 p = double(2),
                 log = double(0, default = 0)) {
    if (length(gamma) != dim(x)[1]) stop("Length of gamma vector does not match length of data.")
    if (dim(p)[1] != dim(x)[1]) stop("Dimension mismatch between x and p matrices.")
    if (dim(p)[2] != dim(x)[2]) stop("Dimension mismatch between x and p matrices.")

    ## x is a year by rep matix
    ProbOccNextTime <- psi1
    ll <- 0
    nyears <- dim(x)[1]
    if(nyears >= 1) {
      for(t in 1:nyears) {
        if(nrep[t] > 0) {
          numObs <- sum(x[t,1:nrep[t]])
          if(numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,1:nrep[t]],
                             size = 1, p = p[t,1:nrep[t]], log = 1)))
          ProbUnoccAndCount <- (1-ProbOccNextTime) * (numObs == 0)
          ProbCount <- ProbOccAndCount + ProbUnoccAndCount
          ProbOccGivenCount <- ProbOccAndCount / ProbCount
          ll <- ll + log(ProbCount)
          if(t < nyears)
              ProbOccNextTime <- ProbOccGivenCount * phi +
                  (1-ProbOccGivenCount) * gamma[t]
        } else {
          ## If there were no observations in year t,
          ## simply propagate probability of occupancy forward
          if (t < nyears)
              ProbOccNextTime <- ProbOccNextTime * phi +
                  (1-ProbOccNextTime) * gamma[t]
        }
      }
    }
    if(log) return(ll)
    else return(exp(ll))
    returnType(double(0))
  }
)

dDynOcc_ss <- nimbleFunction(
  ## DynamicOccupancy removes the z's and muZ's from the model and computes
  ## the probability of all reps over all years for one site, species
  run = function(x = double(2),
                 nrep = double(1),
                 psi1 = double(0),
                 phi = double(0),
                 gamma = double(0),
                 p = double(2),
                 log = double(0, default = 0)) {
    # if (length(gamma) != 1) stop("In dDynOcc_vs gamma must be scalar")
    # if (length(phi) != 1) stop("In dDynOcc_vs phi must be scalar")
    if (dim(p)[1] != dim(x)[1]) stop("Dimension mismatch between x and p matrices.")
    if (dim(p)[2] != dim(x)[2]) stop("Dimension mismatch between x and p matrices.")

    ## x is a year by rep matix
    ProbOccNextTime <- psi1
    ll <- 0
    nyears <- dim(x)[1]
    if(nyears >= 1) {
      for(t in 1:nyears) {
        if(nrep[t] > 0) {
          numObs <- sum(x[t,1:nrep[t]])
          if(numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,1:nrep[t]],
                             size = 1, p = p[t,1:nrep[t]], log = 1)))
          ProbUnoccAndCount <- (1-ProbOccNextTime) * (numObs == 0)
          ProbCount <- ProbOccAndCount + ProbUnoccAndCount
          ProbOccGivenCount <- ProbOccAndCount / ProbCount
          ll <- ll + log(ProbCount)
          if(t < nyears)
              ProbOccNextTime <- ProbOccGivenCount * phi +
                  (1-ProbOccGivenCount) * gamma
        } else {
          ## If there were no observations in year t,
          ## simply propagate probability of occupancy forward
          if (t < nyears)
              ProbOccNextTime <- ProbOccNextTime * phi +
                  (1-ProbOccNextTime) * gamma
        }
      }
    }
    if(log) return(ll)
    else return(exp(ll))
    returnType(double(0))
  }
)


registerDistributions(list(
    dDynOcc_vv = list(
        BUGSdist = "dDynOcc_vv(nrep, psi1, phi, gamma, p)",
        Rdist = "dDynOcc_vv(nrep, psi1, phi, gamma, p)",
        types = c('value = double(2)',
                  'nrep = double(1)',
                  'psi1 = double(0)',
                  'phi = double(1)',
                  'gamma = double(1)',
                  'p = double(2)')),
    dDynOcc_vs = list(
        BUGSdist = "dDynOcc_vs(nrep, psi1, phi, gamma, p)",
        Rdist = "dDynOcc_vs(nrep, psi1, phi, gamma, p)",
        types = c('value = double(2)',
                  'nrep = double(1)',
                  'psi1 = double(0)',
                  'phi = double(1)',
                  'gamma = double(0)',
                  'p = double(2)')),
    dDynOcc_sv = list(
        BUGSdist = "dDynOcc_sv(nrep, psi1, phi, gamma, p)",
        Rdist = "dDynOcc_sv(nrep, psi1, phi, gamma, p)",
        types = c('value = double(2)',
                  'nrep = double(1)',
                  'psi1 = double(0)',
                  'phi = double(0)',
                  'gamma = double(1)',
                  'p = double(2)')),
    dDynOcc_ss = list(
        BUGSdist = "dDynOcc_ss(nrep, psi1, phi, gamma, p)",
        Rdist = "dDynOcc_ss(nrep, psi1, phi, gamma, p)",
        types = c('value = double(2)',
                  'nrep = double(1)',
                  'psi1 = double(0)',
                  'phi = double(0)',
                  'gamma = double(0)',
                  'p = double(2)'))
))
