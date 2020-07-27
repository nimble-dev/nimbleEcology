# nsite <- 50
#
# nobs <- 2 + rpois(nsite, 2)
#
# occ_data <- matrix(data = rnorm(nsite * 4), ncol = 4)
# det_data <- array(dim = c(nsite, max(nobs), 4))
# y <- matrix(ncol = max(nobs), nrow = nsite)
# z <- numeric(nsite)
#
# for (i in 1:nsite) {
#   z[i] <- rbinom(1, 1, nimble::expit(sum(occ_data[i,])))
#   for (j in 1:4) det_data[i, 1:nobs[i], j] <- rnorm(1:nobs[i])
#   y[i, 1:nobs[i]] <- rbinom(nobs[i], z[i], nimble::expit(rowSums(det_data[i,1:nobs[i],])))
# }
# image(y)




#' @name fit_Occ
#' @title Fit an occupancy model with NIMBLE
#' @param y A matrix of dimension (# sites x max # observations at one site).
#'   The element \code{y[i,j]} indicates whether the species was detected (1) or
#'   not detected (0), or that an observation did not take place (NA). Currently
#'   NAs are only allowed at the end of a row, so the observation history [0, 0,
#'   1, NA, NA] is allowed but not [0, 0, NA, 1, NA].
#' @param occ_data A data frame giving the covariate data for the occupancy
#'   process. These are site-level data (shared across visits to a site).
#' @param det_data An array. Each slice \code{det_data[,,i]} gives the
#'   observation-level data for the ith covariate on the detection process, the
#'   dimension of which will match \code{y}.
#'   If only site-level covariates are needed
#' @param method How should the model be estimated? One of "MLE" or "MCMC."
fit_Occ <- function(y, formula_occ, formula_det,
                    occ_data, det_data,
                    method = c("MLE", "MCMC"),
                    return_nimbleModel = FALSE,
                    optim_method = NULL,
                    niter = NULL, nburnin = NULL, nchains = NULL,
                    silent = FALSE) {

  occ_data <- model.matrix(formula_occ, data = data)
  det_data <- model.matrix(formula_det, data = data)

  occ_code <- getNCodeBlock(method = method)

}


# TODO: Flexibility options it could be nice to have
#  - MLE
#    - Choose link fn for occupancy and
#    - Missing obs in the middle of y
#    -
#  - MCMC
#    - Specify informative priors
#    -
#    -


