
# Internal function. Basically a big switch() call to store / retrieve nimbleCode.
getNCodeBlock <- function(model, method, use) {

  if (!(model %in% c("Occ", "Nmixture", "DynOcc", "HMM", "DHMM", "CJS")))
  if (!(method %in% c("MLE", "MCMC"))) stop(paste0("Unknown method ", method))

  if (model == "Occ") {
    if (method == "MLE")
      model_code <- nimbleCode({
        for (o in 1:nobs) {
          logit(p[j]) <- inprod(detDat[o, 1:nDetCovs], detBeta[1:nDetCovs])
        }
        for (i in 1:nsiteRep) {
          logit(psi[i]) <- inprod(occDat[i, 1:nOccCovs], occBeta[1:nOccCovs])
          y[start[i]:end[i]]] ~ dOcc_v(probOcc = psi[i],
                                       probDetect = p[start[i]:end[i]],
                                       len = end[i] - start[i] + 1)
        }
        for (i in (nsiteRep + 1):(nsiteRep + nsiteOneOnly)) {
          logit(psi[i]) <- inprod(occDat[i, 1:nOccCovs], occBeta[1:nOccCovs])
          y[start[i]:end[i]]] ~ dOcc_v(probOcc = psi[i],
                                       probDetect = p[start[i]:end[i]],
                                       len = end[i] - start[i] + 1)
        }
      })

    else (method == "MCMC")
      model_code <- nimbleCode({
        # Include some uninformative priors
      })
  } else {
    stop(paste0("Sorry, I haven't coded ", model, " yet."))
  }


  return(model_code)
}
