# dDynOcc
# To be filled in with dynamic occupancy model distribution(s).
#' Dynamic occupancy distribution for use in NIMBLE models
#'
#' @aliases dDynOcc rOcc rOcc_s dOcc_v rOcc_v
#'
#' @export
dDynamicOccupancy <- nimbleFunction(
  ## DynamicOccupancy removes the z's and muZ's from the model and computes
  ## the probability of all reps over all years for one site, species
  run = function(x = double(2),
                 nrep = double(1),
                 psi1 = double(0),
                 phi = double(1),
                 gamma = double(1),
                 p = double(2),
                 log = double(0, default = 0)) {
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
          if(t < nyears)
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

registerDistributions(list(
    dDynamicOccupancy = list(
        BUGSdist = "dDynamicOccupancy(nrep, psi1, phi, gamma, p)",
        Rdist = "dDynamicOccupancy(nrep, psi1, phi, gamma, p)",
        types = c('value = double(2)',
                  'nrep = double(1)',
                  'psi1 = double(0)',
                  'phi = double(1)',
                  'gamma = double(1)',
                  'p = double(2)'))
))
