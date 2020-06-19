context("Test that functions work when AD is enabled.")

oldDerivOption <- nimbleOptions("experimentalEnableDerivs")
nimbleOptions(experimentalEnableDerivs = TRUE)

nfm <- nimbleFunction(
  setup = function(model, wrt, nodes) {},
  run = function(x = double(1),
                 order = double(1),
                 reset = logical(0, default=FALSE)) {
    values(model, wrt) <<- x
    ans <- nimDerivs(model$calculate(nodes), wrt = wrt,
                     order = order, reset = reset)
    return(ans)
    returnType(ADNimbleList())
  }
)

########################## Test dCJS ############################

test_that("dCJS works with AD", {

# Test code with dCJS_sv
  cjs_code <- nimbleCode({
    logit(pSurv) <- pSurv_int
    for (i in 1:10) {
      logit(pCap[i]) <- pCap_int + beta_pCap * (i - 5.5)

      x[i, 1:ntime] ~ dCJS_sv(probSurvive = pSurv,
                              probCapture = pCap[1:ntime],
                              len = ntime)
    }
  })
  cjs_model <- nimbleModel(cjs_code,
                           constants = list(
                             ntime = 10
                           ), inits = list(
                             beta_pCap = 0.05, pCap_int = 0, pSurv_int = 2
                           ))
  cjs_model$simulate("x")
  cjs_wrt <- cjs_model$expandNodeNames(c("beta_pCap", "pCap_int", "pSurv_int"))
  cjs_nodes <- cjs_model$getDependencies(cjs_wrt)
  cjs_nfm <- nfm(model = cjs_model,
                 wrt = cjs_wrt,
                 nodes = cjs_nodes)
  Ccjs_model <- compileNimble(cjs_model)
  Ccjs_nfm <- compileNimble(cjs_nfm)
  cjs_sv_result <- Ccjs_nfm$run(x = values(cjs_model, cjs_wrt), order = c(0,1,2))
  expect_true(all(!is.na(cjs_sv_result$hessian)))

# Test code with dCJS_vs
  cjs_code <- nimbleCode({
    logit(pCap) <- pCap_int
    for (i in 1:9) logit(pSurv[i]) <- pSurv_int + beta_pSurv * (i - 5.5)
    for (i in 1:10) {

      x[i, 1:ntime] ~ dCJS_vs(probSurvive = pSurv[1:(ntime - 1)],
                              probCapture = pCap,
                              len = ntime)
    }
  })
  cjs_model <- nimbleModel(cjs_code,
                           constants = list(
                             ntime = 10
                           ), inits = list(
                             beta_pSurv = 0.05, pCap_int = 0, pSurv_int = 2
                           ))
  cjs_model$simulate("x")
  cjs_wrt <- cjs_model$expandNodeNames(c("beta_pSurv", "pCap_int", "pSurv_int"))
  cjs_nodes <- cjs_model$getDependencies(cjs_wrt)
  cjs_nfm <- nfm(model = cjs_model,
                 wrt = cjs_wrt,
                 nodes = cjs_nodes)
  Ccjs_model <- compileNimble(cjs_model)
  Ccjs_nfm <- compileNimble(cjs_nfm)
  cjs_sv_result <- Ccjs_nfm$run(x = values(cjs_model, cjs_wrt), order = c(0,1,2))
  expect_true(all(!is.na(cjs_sv_result$hessian)))

# Test code with dCJS_vv
  cjs_code <- nimbleCode({
    for (i in 1:9) logit(pSurv[i]) <- pSurv_int + beta_pSurv * (i - 5.5)
    for (i in 1:10) {
      logit(pCap[i]) <- pCap_int + beta_pCap * (i - 5.5)

      x[i, 1:ntime] ~ dCJS_vv(probSurvive = pSurv[1:(ntime - 1)],
                              probCapture = pCap[1:ntime],
                              len = ntime)
    }
  })
  cjs_model <- nimbleModel(cjs_code,
                           constants = list(
                             ntime = 10
                           ), inits = list(
                             beta_pSurv = 0.05, beta_pCap = 0.1, pCap_int = 0, pSurv_int = 2
                           ))
  cjs_model$simulate("x")
  cjs_wrt <- cjs_model$expandNodeNames(c("beta_pSurv", "pCap_int", "pSurv_int", "beta_pCap"))
  cjs_nodes <- cjs_model$getDependencies(cjs_wrt)
  cjs_nfm <- nfm(model = cjs_model,
                 wrt = cjs_wrt,
                 nodes = cjs_nodes)
  Ccjs_model <- compileNimble(cjs_model)
  Ccjs_nfm <- compileNimble(cjs_nfm)
  cjs_sv_result <- Ccjs_nfm$run(x = values(cjs_model, cjs_wrt), order = c(0,1,2))
  expect_true(all(!is.na(cjs_sv_result$hessian)))

# Test code with dCJS_ss
  cjs_code <- nimbleCode({
    logit(pSurv) <- pSurv_int
    logit(pCap) <- pCap_int
    for (i in 1:10) {
      x[i, 1:ntime] ~ dCJS_ss(probSurvive = pSurv,
                              probCapture = pCap,
                              len = ntime)
    }
  })
  cjs_model <- nimbleModel(cjs_code,
                           constants = list(
                             ntime = 10
                           ), inits = list(
                             pCap_int = 0, pSurv_int = 2
                           ))
  cjs_model$simulate("x")
  cjs_wrt <- cjs_model$expandNodeNames(c("pCap_int", "pSurv_int"))
  cjs_nodes <- cjs_model$getDependencies(cjs_wrt)
  cjs_nfm <- nfm(model = cjs_model,
                 wrt = cjs_wrt,
                 nodes = cjs_nodes)
  Ccjs_model <- compileNimble(cjs_model)
  Ccjs_nfm <- compileNimble(cjs_nfm)
  cjs_sv_result <- Ccjs_nfm$run(x = values(cjs_model, cjs_wrt), order = c(0,1,2))
  expect_true(all(!is.na(cjs_sv_result$hessian)))
})

########################## Test dOcc ############################
test_that("dOcc works with AD", {

# Test dOcc_v
  occ_code <- nimbleCode({
    for (i in 1:nsite) {
      logit(psi[i]) <- inprod(psi_beta[1:3], occu_cov[i, 1:3])
      for (j in 1:nvisit) {
        logit(p[i,j]) <- inprod(p_beta[1:3], detect_cov[i, j, 1:3])
      }
      y[i, 1:nvisit] ~ dOcc_v(probOcc = psi[i],
                              probDetect = p[i, 1:nvisit],
                              len = nvisit)
    }
  })
  nsite <- 30
  nvisit <- 3
  detect_cov <- array(rnorm(nsite * nvisit * 3),
                      dim = c(nsite, nvisit, 3))
  detect_cov[,,1] <- 1
  occu_cov <- matrix(data = rnorm(nsite*3), nrow = nsite)
  occu_cov[,1] <- 1
  psi_beta <- c(0, 1, -1)
  p_beta <- c(1, 1, -1)
  occ_model <- nimbleModel(code = occ_code,
                           constants = list(
                             nsite = nsite,
                             nvisit = nvisit),
                           data = list(
                             occu_cov = occu_cov,
                             detect_cov = detect_cov
                           ),
                           inits = list(
                             psi_beta = psi_beta,
                             p_beta = p_beta
                           ))
  occ_model$simulate("y")
  C_occ_model <- compileNimble(occ_model)
  wrt <- c(occ_model$expandNodeNames("psi_beta"),
           occ_model$expandNodeNames("p_beta"))
  nodes <- occ_model$getDependencies(wrt)
  nfm1 <- nfm(occ_model, wrt, nodes)
  Cnfm1 <- compileNimble(nfm1)
  occ_result <- Cnfm1$run(x = rep(0, 6), order = c(0,1,2))
  expect_true(all(!is.na(occ_result$hessian)))

# Test dOcc_s
  occ_code <- nimbleCode({
    for (i in 1:nsite) {
      logit(psi[i]) <- inprod(psi_beta[1:3], occu_cov[i, 1:3])
      logit(p[i]) <- inprod(p_beta[1:3], detect_cov[i, 1:3])
      y[i, 1:nvisit] ~ dOcc_s(probOcc = psi[i],
                              probDetect = p[i],
                              len = nvisit)
    }
  })
  nsite <- 30
  nvisit <- 3
  detect_cov <- matrix(data = rnorm(nsite*3), nrow = nsite)
  detect_cov[,1] <- 1
  occu_cov <- matrix(data = rnorm(nsite*3), nrow = nsite)
  occu_cov[,1] <- 1
  psi_beta <- c(0, 1, -1)
  p_beta <- c(1, 1, -1)
  occ_model <- nimbleModel(code = occ_code,
                           constants = list(
                             nsite = nsite,
                             nvisit = nvisit),
                           data = list(
                             occu_cov = occu_cov,
                             detect_cov = detect_cov
                           ),
                           inits = list(
                             psi_beta = psi_beta,
                             p_beta = p_beta
                           ))
  occ_model$simulate("y")
  C_occ_model <- compileNimble(occ_model)
  wrt <- c(occ_model$expandNodeNames("psi_beta"),
           occ_model$expandNodeNames("p_beta"))
  nodes <- occ_model$getDependencies(wrt)
  nfm1 <- nfm(occ_model, wrt, nodes)
  Cnfm1 <- compileNimble(nfm1)
  occ_result <- Cnfm1$run(x = rep(0, 6), order = c(0,1,2))
  expect_true(all(!is.na(occ_result$hessian)))
})


########################## Test dHMM ############################
test_that("dHMM works with AD", {

# Test dHMM
  hmm_code <- nimbleCode({
    for (i in 1:10) {
      x[i, 1:ntime] ~ dHMM(init = inits[1:nstate],
                           probObs = pO[1:nstate, 1:nobs],
                           probTrans = pT[1:nstate, 1:nstate],
                           len = ntime,
                           checkRowSums = 1)
    }
  })
  hmm_model <- nimbleModel(hmm_code,
                           constants = list(
                             ntime = 10,
                             nstate = 3,
                             nobs = 2
                           ), inits = list(
                             inits = c(0.9, 0.1, 0),
                             pO = matrix(c(0.9, 0.1,
                                           0.8, 0.2,
                                           0, 1), nrow = 3, byrow = TRUE),
                             pT = matrix(c(0.8, 0.2, 0,
                                           0, 0.7, 0.3,
                                           0, 0, 1), nrow = 3, byrow = TRUE)))
  hmm_model$simulate("x")
  hmm_model$x
  hmm_wrt <- hmm_model$expandNodeNames(c("inits", "pO", "pT"))
  hmm_nodes <- hmm_model$getDependencies(hmm_wrt)
  hmm_nfm <- nfm(model = hmm_model,
                 wrt = hmm_wrt,
                 nodes = hmm_nodes)
  Chmm_model <- compileNimble(hmm_model)
  Chmm_nfm <- compileNimble(hmm_nfm)
  hmm_result <- Chmm_nfm$run(x = values(hmm_model, hmm_wrt), order = c(0,1,2))
  expect_true(all(!is.na(hmm_result$hessian)))

  # Test dHMMo
  hmm_code <- nimbleCode({
    for (i in 1:10) {
      x[i, 1:ntime] ~ dHMMo(init = inits[1:nstate],
                            probObs = pO[1:nstate, 1:nobs, 1:ntime],
                            probTrans = pT[1:nstate, 1:nstate],
                            len = ntime,
                            checkRowSums = 1)
    }
  })

  hmm_model <- nimbleModel(hmm_code,
                           constants = list(
                             ntime = 10,
                             nstate = 3,
                             nobs = 2
                           ), inits = list(
                             inits = c(0.9, 0.1, 0),
                             pO = array(rep(c(0.9, 0.8, 0.1,
                                              0.1, 0.2, 0.9), 10), dim = c(3, 2, 10)),
                             pT = matrix(c(0.8, 0.2, 0,
                                           0, 0.7, 0.3,
                                           0, 0, 1), nrow = 3, byrow = TRUE)))
  hmm_model$simulate("x")
  hmm_model$x
  hmm_model$calculate()
  Chmm_model <- compileNimble(hmm_model)
  hmm_wrt <- hmm_model$expandNodeNames(c("inits", "pO", "pT"))
  hmm_nodes <- hmm_model$getDependencies(hmm_wrt)
  hmm_nfm <- nfm(model = hmm_model,
                 wrt = hmm_wrt,
                 nodes = hmm_nodes)

  Chmm_nfm <- compileNimble(hmm_nfm)
  hmm_result <- Chmm_nfm$run(x = values(hmm_model, hmm_wrt), order = c(0,1,2))
  expect_true(all(!is.na(hmm_result$hessian)))
})



############################## Test dDHMM #####################################
test_that("dHMM works with AD", {

  # Test DHMM
  dhmm_code <- nimbleCode({
    for (i in 1:10) {
      x[i, 1:ntime] ~ dDHMM(init = inits[1:nstate],
                            probObs = pO[1:nstate, 1:nobs],
                            probTrans = pT[1:nstate, 1:nstate, 1:(ntime-1)],
                            len = ntime,
                            checkRowSums = 1)
    }
  })
  dhmm_model <- nimbleModel(dhmm_code,
                            constants = list(
                              ntime = 10,
                              nstate = 3,
                              nobs = 2
                            ), inits = list(
                              inits = c(0.9, 0.1, 0),
                              pO = matrix(c(0.9, 0.1,
                                            0.8, 0.2,
                                            0, 1), nrow = 3, byrow = TRUE),
                              pT = array(rep(c(0.8, 0, 0,
                                               0.2, 0.7, 0,
                                               0, 0.3, 1), 9), dim = c(3, 3, 9))))
  dhmm_model$simulate("x")
  dhmm_wrt <- dhmm_model$expandNodeNames(c("inits", "pO", "pT"))
  dhmm_nodes <- dhmm_model$getDependencies(dhmm_wrt)
  dhmm_nfm <- nfm(model = dhmm_model,
                  wrt = dhmm_wrt,
                  nodes = dhmm_nodes)
  Cdhmm_model <- compileNimble(dhmm_model)
  Cdhmm_nfm <- compileNimble(dhmm_nfm)
  dhmm_result <- Cdhmm_nfm$run(x = values(dhmm_model, dhmm_wrt), order = c(0,1,2))
  expect_true(all(!is.na(dhmm_result$hessian)))

  # Test DHMMo
  # Test DHMM
  dhmm_code <- nimbleCode({
    for (i in 1:10) {
      x[i, 1:ntime] ~ dDHMMo(init = inits[1:nstate],
                            probObs = pO[1:nstate, 1:nobs, 1:ntime],
                            probTrans = pT[1:nstate, 1:nstate, 1:(ntime-1)],
                            len = ntime,
                            checkRowSums = 1)
    }
  })
  dhmm_model <- nimbleModel(dhmm_code,
                            constants = list(
                              ntime = 10,
                              nstate = 3,
                              nobs = 2
                            ), inits = list(
                              inits = c(0.9, 0.1, 0),
                              pO = array(rep(c(0.9, 0.8, 0.1,
                                               0.1, 0.2, 0.9), 10), dim = c(3, 2, 10)),
                              pT = array(rep(c(0.8, 0, 0,
                                               0.2, 0.7, 0,
                                               0, 0.3, 1), 9), dim = c(3, 3, 9))))
  dhmm_model$simulate("x")
  dhmm_wrt <- dhmm_model$expandNodeNames(c("inits", "pO", "pT"))
  dhmm_nodes <- dhmm_model$getDependencies(dhmm_wrt)
  dhmm_nfm <- nfm(model = dhmm_model,
                  wrt = dhmm_wrt,
                  nodes = dhmm_nodes)
  Cdhmm_model <- compileNimble(dhmm_model)
  Cdhmm_nfm <- compileNimble(dhmm_nfm)
  dhmm_result <- Cdhmm_nfm$run(x = values(dhmm_model, dhmm_wrt), order = c(0,1,2))
  expect_true(all(!is.na(dhmm_result$hessian)))

})


########################### Test dDynOcc_ss* ##################################
# Since there are so many flavors of dDynOcc, I'm going to split them up
# across multiple test_that blocks.
test_that("dDynOcc_ss* works with AD", {

# Test DynOcc_sss
  dynocc_code <- nimbleCode({
    x[1:nssn, 1:nvisit] ~ dDynOcc_sss(inits,
                                     probCol,
                                     probPer,
                                     p = p,
                                     start = start[1:nssn],
                                     end = end[1:nssn])

    logit(p) <- p_int
    logit(probCol) <- col_int
    logit(probPer) <- per_int
  })
  dynocc_model <- nimbleModel(code = dynocc_code,
                              constants = list(
                                nvisit = 4,
                                nssn = 3,
                                start = rep(1, 10),
                                end = rep(4, 10)
                              ), inits = list(
                                p_int = 0.5,
                                col_int = 0.5,
                                per_int = -0.5,
                                inits = 0.9
                              ))
  dynocc_model$simulate("x")
  dynocc_model$x
  dynocc_model$calculate()
  dynocc_nfm <- nfm(model = dynocc_model,
                    wrt = c("p_int", "col_int", "per_int", "inits"),
                    nodes = dynocc_model$getDependencies(c("p_int", "col_int",
                                                           "per_int", "inits")))
  Cdynocc_model <- compileNimble(dynocc_model)
  Cdynocc_nfm <- compileNimble(dynocc_nfm)
  dynocc_result <- Cdynocc_nfm$run(x = rep(0.5, 4), order = c(0,1,2))
  expect_true(all(!is.na(dynocc_result$hessian)))

# Test DynOcc_ssv
  dynocc_code <- nimbleCode({
    x[1:nssn, 1:nvisit] ~ dDynOcc_ssv(inits,
                                     probCol,
                                     probPer,
                                     p = p[1:nssn],
                                     start = start[1:nssn],
                                     end = end[1:nssn])

    for (i in 1:nssn) logit(p[i]) <- p_int + i - 2
    logit(probCol) <- col_int
    logit(probPer) <- per_int
  })
  dynocc_model <- nimbleModel(code = dynocc_code,
                              constants = list(
                                nvisit = 4,
                                nssn = 3,
                                start = rep(1, 3),
                                end = rep(4, 3)
                              ), inits = list(
                                p_int = 0.5,
                                col_int = 0.5,
                                per_int = -0.5,
                                inits = 0.9
                              ))
  dynocc_model$simulate("x[1:3, 1:4]", includeData = TRUE)
  dynocc_model$x
  dynocc_model$calculate()
  dynocc_nfm <- nfm(model = dynocc_model,
                    wrt = c("p_int", "col_int", "per_int", "inits"),
                    nodes = dynocc_model$getDependencies(c("p_int", "col_int",
                                                           "per_int", "inits")))
  Cdynocc_model <- compileNimble(dynocc_model)
  Cdynocc_nfm <- compileNimble(dynocc_nfm)
  dynocc_result <- Cdynocc_nfm$run(x = rep(0.5, 4), order = c(0,1,2))
  expect_true(all(!is.na(dynocc_result$hessian)))

# Test DynOcc_ssm
  dynocc_code <- nimbleCode({
    x[1:nssn, 1:nvisit] ~ dDynOcc_ssm(inits,
                                      probCol,
                                      probPer,
                                      p = p[1:nssn, 1:nvisit],
                                      start = start[1:nssn],
                                      end = end[1:nssn])

    for (i in 1:nssn){
      for (j in 1:nvisit) {
        logit(p[i, j]) <- p_int + i - 2 + j - 2
      }
    }
    logit(probCol) <- col_int
    logit(probPer) <- per_int
  })
  dynocc_model <- nimbleModel(code = dynocc_code,
                              constants = list(
                                nvisit = 4,
                                nssn = 3,
                                start = rep(1, 3),
                                end = rep(4, 3)
                              ), inits = list(
                                p_int = 0.5,
                                col_int = 0.5,
                                per_int = -0.5,
                                inits = 0.9
                              ))
  dynocc_model$simulate("x[1:3, 1:4]", includeData = TRUE)
  dynocc_model$x
  dynocc_model$calculate()
  dynocc_nfm <- nfm(model = dynocc_model,
                    wrt = c("p_int", "col_int", "per_int", "inits"),
                    nodes = dynocc_model$getDependencies(c("p_int", "col_int",
                                                           "per_int", "inits")))
  Cdynocc_model <- compileNimble(dynocc_model)
  Cdynocc_nfm <- compileNimble(dynocc_nfm)
  dynocc_result <- Cdynocc_nfm$run(x = rep(0.5, 4), order = c(0,1,2))
  expect_true(all(!is.na(dynocc_result$hessian)))

})

########################### Test dDynOcc_sv* ##################################
# Since there are so many flavors of dDynOcc, I'm going to split them up
# across multiple test_that blocks.
test_that("dDynOcc_sv* works with AD", {

  # Test DynOcc_svs
  dynocc_code <- nimbleCode({
    x[1:nssn, 1:nvisit] ~ dDynOcc_svs(inits,
                                     probCol,
                                     probPer[1:(nssn-1)],
                                     p = p,
                                     start = start[1:nssn],
                                     end = end[1:nssn])

    logit(p) <- p_int
    logit(probCol) <- col_int
    for (i in 1:(nssn - 1)) logit(probPer[i]) <- per_int + i - 1
  })
  dynocc_model <- nimbleModel(code = dynocc_code,
                              constants = list(
                                nvisit = 4,
                                nssn = 3,
                                start = rep(1, 10),
                                end = rep(4, 10)
                              ), inits = list(
                                p_int = 0.5,
                                col_int = 0.5,
                                per_int = -0.5,
                                inits = 0.9
                              ))
  dynocc_model$simulate("x")
  dynocc_model$x
  dynocc_model$calculate()
  dynocc_nfm <- nfm(model = dynocc_model,
                    wrt = c("p_int", "col_int", "per_int", "inits"),
                    nodes = dynocc_model$getDependencies(c("p_int", "col_int",
                                                           "per_int", "inits")))
  Cdynocc_model <- compileNimble(dynocc_model)
  Cdynocc_nfm <- compileNimble(dynocc_nfm)
  dynocc_result <- Cdynocc_nfm$run(x = rep(0.5, 4), order = c(0,1,2))
  expect_true(all(!is.na(dynocc_result$hessian)))

  # Test DynOcc_svv
  dynocc_code <- nimbleCode({
    x[1:nssn, 1:nvisit] ~ dDynOcc_svv(inits,
                                      probCol,
                                      probPer[1:(nssn-1)],
                                      p = p[1:nssn],
                                      start = start[1:nssn],
                                      end = end[1:nssn])

    for (i in 1:nssn) logit(p[i]) <- p_int + i - 2
    logit(probCol) <- col_int
    for (i in 1:(nssn - 1)) logit(probPer[i]) <- per_int + i - 1
  })
  dynocc_model <- nimbleModel(code = dynocc_code,
                              constants = list(
                                nvisit = 4,
                                nssn = 3,
                                start = rep(1, 3),
                                end = rep(4, 3)
                              ), inits = list(
                                p_int = 0.5,
                                col_int = 0.5,
                                per_int = -0.5,
                                inits = 0.9
                              ))
  dynocc_model$simulate("x[1:3, 1:4]", includeData = TRUE)
  dynocc_model$x
  dynocc_model$calculate()
  dynocc_nfm <- nfm(model = dynocc_model,
                    wrt = c("p_int", "col_int", "per_int", "inits"),
                    nodes = dynocc_model$getDependencies(c("p_int", "col_int",
                                                           "per_int", "inits")))
  Cdynocc_model <- compileNimble(dynocc_model)
  Cdynocc_nfm <- compileNimble(dynocc_nfm)
  dynocc_result <- Cdynocc_nfm$run(x = rep(0.5, 4), order = c(0,1,2))
  expect_true(all(!is.na(dynocc_result$hessian)))

  # Test DynOcc_svm
  dynocc_code <- nimbleCode({
    x[1:nssn, 1:nvisit] ~ dDynOcc_svm(inits,
                                      probCol,
                                      probPer[1:(nssn-1)],
                                      p = p[1:nssn, 1:nvisit],
                                      start = start[1:nssn],
                                      end = end[1:nssn])

    for (i in 1:nssn){
      for (j in 1:nvisit) {
        logit(p[i, j]) <- p_int + i - 2 + j - 2
      }
    }
    logit(probCol) <- col_int
    for (i in 1:(nssn - 1)) logit(probPer[i]) <- per_int + i - 1
  })
  dynocc_model <- nimbleModel(code = dynocc_code,
                              constants = list(
                                nvisit = 4,
                                nssn = 3,
                                start = rep(1, 3),
                                end = rep(4, 3)
                              ), inits = list(
                                p_int = 0.5,
                                col_int = 0.5,
                                per_int = -0.5,
                                inits = 0.9
                              ))
  dynocc_model$simulate("x[1:3, 1:4]", includeData = TRUE)
  dynocc_model$x
  dynocc_model$calculate()
  dynocc_nfm <- nfm(model = dynocc_model,
                    wrt = c("p_int", "col_int", "per_int", "inits"),
                    nodes = dynocc_model$getDependencies(c("p_int", "col_int",
                                                           "per_int", "inits")))
  Cdynocc_model <- compileNimble(dynocc_model)
  Cdynocc_nfm <- compileNimble(dynocc_nfm)
  dynocc_result <- Cdynocc_nfm$run(x = rep(0.5, 4), order = c(0,1,2))
  expect_true(all(!is.na(dynocc_result$hessian)))

})

########################### Test dDynOcc_vv* ##################################
# Since there are so many flavors of dDynOcc, I'm going to split them up
# across multiple test_that blocks.
test_that("dDynOcc_vv* works with AD", {

  # Test DynOcc_vvs
  dynocc_code <- nimbleCode({
    x[1:nssn, 1:nvisit] ~ dDynOcc_vvs(inits,
                                      probCol[1:(nssn-1)],
                                      probPer[1:(nssn-1)],
                                      p = p,
                                      start = start[1:nssn],
                                      end = end[1:nssn])

    logit(p) <- p_int
    for (i in 1:(nssn - 1)) logit(probCol[i]) <- col_int + i - 1
    for (i in 1:(nssn - 1)) logit(probPer[i]) <- per_int + i - 1
  })
  dynocc_model <- nimbleModel(code = dynocc_code,
                              constants = list(
                                nvisit = 4,
                                nssn = 3,
                                start = rep(1, 10),
                                end = rep(4, 10)
                              ), inits = list(
                                p_int = 0.5,
                                col_int = 0.5,
                                per_int = -0.5,
                                inits = 0.9
                              ))
  dynocc_model$simulate("x")
  dynocc_model$x
  dynocc_model$calculate()
  dynocc_nfm <- nfm(model = dynocc_model,
                    wrt = c("p_int", "col_int", "per_int", "inits"),
                    nodes = dynocc_model$getDependencies(c("p_int", "col_int",
                                                           "per_int", "inits")))
  Cdynocc_model <- compileNimble(dynocc_model)
  Cdynocc_nfm <- compileNimble(dynocc_nfm)
  dynocc_result <- Cdynocc_nfm$run(x = rep(0.5, 4), order = c(0,1,2))
  expect_true(all(!is.na(dynocc_result$hessian)))

  # Test DynOcc_vvv
  dynocc_code <- nimbleCode({
    x[1:nssn, 1:nvisit] ~ dDynOcc_vvv(inits,
                                      probCol[1:(nssn-1)],
                                      probPer[1:(nssn-1)],
                                      p = p[1:nssn],
                                      start = start[1:nssn],
                                      end = end[1:nssn])

    for (i in 1:nssn) logit(p[i]) <- p_int + i - 2
    for (i in 1:(nssn - 1)) logit(probCol[i]) <- col_int + i - 1
    for (i in 1:(nssn - 1)) logit(probPer[i]) <- per_int + i - 1
  })
  dynocc_model <- nimbleModel(code = dynocc_code,
                              constants = list(
                                nvisit = 4,
                                nssn = 3,
                                start = rep(1, 3),
                                end = rep(4, 3)
                              ), inits = list(
                                p_int = 0.5,
                                col_int = 0.5,
                                per_int = -0.5,
                                inits = 0.9
                              ))
  dynocc_model$simulate("x[1:3, 1:4]", includeData = TRUE)
  dynocc_model$x
  dynocc_model$calculate()
  dynocc_nfm <- nfm(model = dynocc_model,
                    wrt = c("p_int", "col_int", "per_int", "inits"),
                    nodes = dynocc_model$getDependencies(c("p_int", "col_int",
                                                           "per_int", "inits")))
  Cdynocc_model <- compileNimble(dynocc_model)
  Cdynocc_nfm <- compileNimble(dynocc_nfm)
  dynocc_result <- Cdynocc_nfm$run(x = rep(0.5, 4), order = c(0,1,2))
  expect_true(all(!is.na(dynocc_result$hessian)))

  # Test DynOcc_vvm
  dynocc_code <- nimbleCode({
    x[1:nssn, 1:nvisit] ~ dDynOcc_vvm(inits,
                                      probCol[1:(nssn-1)],
                                      probPer[1:(nssn-1)],
                                      p = p[1:nssn, 1:nvisit],
                                      start = start[1:nssn],
                                      end = end[1:nssn])

    for (i in 1:nssn){
      for (j in 1:nvisit) {
        logit(p[i, j]) <- p_int + i - 2 + j - 2
      }
    }
    for (i in 1:(nssn - 1)) logit(probCol[i]) <- col_int + i - 1
    for (i in 1:(nssn - 1)) logit(probPer[i]) <- per_int + i - 1
  })
  dynocc_model <- nimbleModel(code = dynocc_code,
                              constants = list(
                                nvisit = 4,
                                nssn = 3,
                                start = rep(1, 3),
                                end = rep(4, 3)
                              ), inits = list(
                                p_int = 0.5,
                                col_int = 0.5,
                                per_int = -0.5,
                                inits = 0.9
                              ))
  dynocc_model$simulate("x[1:3, 1:4]", includeData = TRUE)
  dynocc_model$x
  dynocc_model$calculate()
  dynocc_nfm <- nfm(model = dynocc_model,
                    wrt = c("p_int", "col_int", "per_int", "inits"),
                    nodes = dynocc_model$getDependencies(c("p_int", "col_int",
                                                           "per_int", "inits")))
  Cdynocc_model <- compileNimble(dynocc_model)
  Cdynocc_nfm <- compileNimble(dynocc_nfm)
  dynocc_result <- Cdynocc_nfm$run(x = rep(0.5, 4), order = c(0,1,2))
  expect_true(all(!is.na(dynocc_result$hessian)))

})


########################### Test dDynOcc_vs* ##################################
# Since there are so many flavors of dDynOcc, I'm going to split them up
# across multiple test_that blocks.
test_that("dDynOcc_vs* works with AD", {

  # Test DynOcc_vvs
  dynocc_code <- nimbleCode({
    x[1:nssn, 1:nvisit] ~ dDynOcc_vss(inits,
                                      probCol[1:(nssn-1)],
                                      probPer,
                                      p = p,
                                      start = start[1:nssn],
                                      end = end[1:nssn])

    logit(p) <- p_int
    for (i in 1:(nssn - 1)) logit(probCol[i]) <- col_int + i - 1
    logit(probPer) <- per_int
  })
  dynocc_model <- nimbleModel(code = dynocc_code,
                              constants = list(
                                nvisit = 4,
                                nssn = 3,
                                start = rep(1, 10),
                                end = rep(4, 10)
                              ), inits = list(
                                p_int = 0.5,
                                col_int = 0.5,
                                per_int = -0.5,
                                inits = 0.9
                              ))
  dynocc_model$simulate("x")
  dynocc_model$x
  dynocc_model$calculate()
  dynocc_nfm <- nfm(model = dynocc_model,
                    wrt = c("p_int", "col_int", "per_int", "inits"),
                    nodes = dynocc_model$getDependencies(c("p_int", "col_int",
                                                           "per_int", "inits")))
  Cdynocc_model <- compileNimble(dynocc_model)
  Cdynocc_nfm <- compileNimble(dynocc_nfm)
  dynocc_result <- Cdynocc_nfm$run(x = rep(0.5, 4), order = c(0,1,2))
  expect_true(all(!is.na(dynocc_result$hessian)))

  # Test DynOcc_vsv
  dynocc_code <- nimbleCode({
    x[1:nssn, 1:nvisit] ~ dDynOcc_vsv(inits,
                                      probCol[1:(nssn-1)],
                                      probPer,
                                      p = p[1:nssn],
                                      start = start[1:nssn],
                                      end = end[1:nssn])

    for (i in 1:nssn) logit(p[i]) <- p_int + i - 2
    for (i in 1:(nssn - 1)) logit(probCol[i]) <- col_int + i - 1
    logit(probPer) <- per_int
  })
  dynocc_model <- nimbleModel(code = dynocc_code,
                              constants = list(
                                nvisit = 4,
                                nssn = 3,
                                start = rep(1, 3),
                                end = rep(4, 3)
                              ), inits = list(
                                p_int = 0.5,
                                col_int = 0.5,
                                per_int = -0.5,
                                inits = 0.9
                              ))
  dynocc_model$simulate("x[1:3, 1:4]", includeData = TRUE)
  dynocc_model$x
  dynocc_model$calculate()
  dynocc_nfm <- nfm(model = dynocc_model,
                    wrt = c("p_int", "col_int", "per_int", "inits"),
                    nodes = dynocc_model$getDependencies(c("p_int", "col_int",
                                                           "per_int", "inits")))
  Cdynocc_model <- compileNimble(dynocc_model)
  Cdynocc_nfm <- compileNimble(dynocc_nfm)
  dynocc_result <- Cdynocc_nfm$run(x = rep(0.5, 4), order = c(0,1,2))
  expect_true(all(!is.na(dynocc_result$hessian)))

  # Test DynOcc_vsm
  dynocc_code <- nimbleCode({
    x[1:nssn, 1:nvisit] ~ dDynOcc_vsm(inits,
                                      probCol[1:(nssn-1)],
                                      probPer,
                                      p = p[1:nssn, 1:nvisit],
                                      start = start[1:nssn],
                                      end = end[1:nssn])

    for (i in 1:nssn){
      for (j in 1:nvisit) {
        logit(p[i, j]) <- p_int + i - 2 + j - 2
      }
    }
    for (i in 1:(nssn - 1)) logit(probCol[i]) <- col_int + i - 1
    logit(probPer) <- per_int
  })
  dynocc_model <- nimbleModel(code = dynocc_code,
                              constants = list(
                                nvisit = 4,
                                nssn = 3,
                                start = rep(1, 3),
                                end = rep(4, 3)
                              ), inits = list(
                                p_int = 0.5,
                                col_int = 0.5,
                                per_int = -0.5,
                                inits = 0.9
                              ))
  dynocc_model$simulate("x[1:3, 1:4]", includeData = TRUE)
  dynocc_model$x
  dynocc_model$calculate()
  dynocc_nfm <- nfm(model = dynocc_model,
                    wrt = c("p_int", "col_int", "per_int", "inits"),
                    nodes = dynocc_model$getDependencies(c("p_int", "col_int",
                                                           "per_int", "inits")))
  Cdynocc_model <- compileNimble(dynocc_model)
  Cdynocc_nfm <- compileNimble(dynocc_nfm)
  dynocc_result <- Cdynocc_nfm$run(x = rep(0.5, 4), order = c(0,1,2))
  expect_true(all(!is.na(dynocc_result$hessian)))

})



########### dNmixture not yet implemented with AD ############
# END
