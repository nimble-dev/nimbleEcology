# Test the Occupancy distribution nimbleFunction.

# -----------------------------------------------------------------------------
# 0. Load
# Test scalar-scalar version
test_that("dOcc_s and rOcc_s work", {
  x <- c(1,0,1,1,0)
  probOcc <- 0.4
  probDetect <- 0.7

  probX <- dOcc_s(x, probOcc, probDetect)
  correctProbX <-
    probOcc * probDetect^3 * (1 - probDetect)^2
  expect_equal(probX, correctProbX)

  x2 <- c(0, 0, 0, 0, 0)
  probX2 <- dOcc_s(x2, probOcc, probDetect)
  correctProbX2 <-
    probOcc * (1 - probDetect)^5 + (1-probOcc)
  expect_equal(probX2, correctProbX2)

  lProbX <- dOcc_s(x, probOcc, probDetect, log = TRUE)
  lCorrectProbX <- log(correctProbX)
  expect_equal(lProbX, lCorrectProbX)

  # we must wrap the call to avoid the error that
  # dOcc_s has no setup code but buildDerivs=TRUE so can't be compiled alone
  call_dOcc_s <- nimbleFunction(
    run = function(x = double(1), probOcc=double(),
                   probDetect=double(),
                   log = integer(0, default = 0)) {
      return(dOcc_s(x, probOcc, probDetect, 0, log))
      returnType(double())
    }
  )
  CdOcc_s <- compileNimble(call_dOcc_s)
  CprobX <- CdOcc_s(x, probOcc, probDetect)
  expect_equal(CprobX, probX)

  ClProbX <- CdOcc_s(x, probOcc, probDetect, log = TRUE)
  expect_equal(ClProbX, lProbX)

  set.seed(1)
  nSim <- 10
  xSim <- matrix(nrow = nSim, ncol = 5)
  for(i in 1:nSim)
    xSim[i,] <- rOcc_s(1, probOcc, probDetect, len = 5)
  set.seed(1)
  CrOcc_s <- compileNimble(rOcc_s)
  CxSim <- matrix(nrow = nSim, ncol = 5)
  for(i in 1:nSim)
    CxSim[i,] <- CrOcc_s(1, probOcc, probDetect, len = 5)
  expect_identical(xSim, CxSim)

  nc <- nimbleCode({
    x[1:5] ~ dOcc_s(probOcc, probDetect, len = 5)
    probDetect ~ dunif(0,1)
    probOcc ~ dunif(0,1)
  })
  m <- nimbleModel(nc, data = list(x = x),
                   inits = list(probOcc = probOcc,
                                probDetect = probDetect))
  m$calculate()
  MlProbX <- m$getLogProb("x")
  expect_equal(MlProbX, lProbX)

  cm <- compileNimble(m)
  cm$calculate()
  CMlProbX <- cm$getLogProb("x")
  expect_equal(CMlProbX, lProbX)

  simNodes <- m$getDependencies(c('probOcc', 'probDetect'), self = FALSE)
  mxSim <- matrix(nrow = nSim, ncol = 5)
  set.seed(1)
  for(i in 1:nSim) {
    m$simulate(simNodes, includeData = TRUE)
    mxSim[i,] <- m$x
  }
  expect_identical(mxSim, xSim)

  CmxSim <- matrix(nrow = nSim, ncol = 5)
  set.seed(1)
  for(i in 1:nSim) {
    cm$simulate(simNodes, includeData = TRUE)
    CmxSim[i,] <- cm$x
  }
  expect_identical(CmxSim, mxSim)

# Test imputing value for all NAs
  xNA <- rep(NA, length(x))
  mNA <- nimbleModel(nc, data = list(x = xNA),
         inits = list(probOcc = probOcc,
                      probDetect = probDetect))
  mNAConf <- configureMCMC(mNA)
  mNAConf$addMonitors('x')
  mNA_MCMC <- buildMCMC(mNAConf)
  cmNA <- compileNimble(mNA, mNA_MCMC)
  set.seed(0)
  cmNA$mNA_MCMC$run(10)
# Did the imputed values come back?
  expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x[1]"])))
})


# Test v version
test_that("dOcc_v works", {
  x <- c(1,0,1,1,0)
  probOcc <- 0.4
  probDetect <- c(0.7, 0.3, 0.5, 0.7, 0.25)

  probX <- dOcc_v(x, probOcc, probDetect)
  correctProbX <-
    probOcc * prod(probDetect[x == 1]) *
    prod(1 - probDetect[x != 1])
  expect_equal(probX, correctProbX)

  x2 <- c(0, 0, 0, 0, 0)
  probX2 <- dOcc_v(x2, probOcc, probDetect)
  correctProbX2 <-
    probOcc * prod(probDetect[x2 == 1]) *
    prod(1 - probDetect[x2 != 1]) + (1-probOcc)
  expect_equal(probX2, correctProbX2)

  lProbX <- dOcc_v(x, probOcc, probDetect, log = TRUE)
  lCorrectProbX <- log(correctProbX)
  expect_equal(lProbX, lCorrectProbX)

  # we must wrap the call to avoid the error that
  # dOcc_v has no setup code but buildDerivs=TRUE so can't be compiled alone
  call_dOcc_v <- nimbleFunction(
    run = function(x = double(1), probOcc=double(0),
                   probDetect=double(1),
                   log = integer(0, default = 0)) {
      return(dOcc_v(x, probOcc, probDetect, 0, log))
      returnType(double())
    }
  )

  CdOcc_v <- compileNimble(call_dOcc_v)
  CprobX <- CdOcc_v(x, probOcc, probDetect)
  expect_equal(CprobX, probX)

  ClProbX <- CdOcc_v(x, probOcc, probDetect, log = TRUE)
  expect_equal(ClProbX, lProbX)

  set.seed(1)
  nSim <- 10
  xSim <- matrix(nrow = nSim, ncol = 5)
  for(i in 1:nSim)
    xSim[i,] <- rOcc_v(1, probOcc, probDetect, len = 5)
  set.seed(1)
  CrOcc_v <- compileNimble(rOcc_v)
  CxSim <- matrix(nrow = nSim, ncol = 5)
  for(i in 1:nSim)
    CxSim[i,] <- CrOcc_v(1, probOcc, probDetect, len = 5)
  expect_identical(xSim, CxSim)

  nc <- nimbleCode({
    x[1:5] ~ dOcc_v(probOcc, probDetect[1:5], len = 5)
    for (i in 1:5) {
      probDetect[i] ~ dunif(0,1)
    }
    probOcc ~ dunif(0,1)
  })
  m <- nimbleModel(nc, data = list(x = x),
                   inits = list(probOcc = probOcc,
                                probDetect = probDetect))
  m$calculate()
  MlProbX <- m$getLogProb("x")
  expect_equal(MlProbX, lProbX)

  cm <- compileNimble(m)
  cm$calculate()
  CMlProbX <- cm$getLogProb("x")
  expect_equal(CMlProbX, lProbX)

  simNodes <- m$getDependencies(c('probOcc', 'probDetect'), self = FALSE)
  mxSim <- matrix(nrow = nSim, ncol = 5)
  set.seed(1)
  for(i in 1:nSim) {
    m$simulate(simNodes, includeData = TRUE)
    mxSim[i,] <- m$x
  }
  expect_identical(mxSim, xSim)

  CmxSim <- matrix(nrow = nSim, ncol = 5)
  set.seed(1)
  for(i in 1:nSim) {
    cm$simulate(simNodes, includeData = TRUE)
    CmxSim[i,] <- cm$x
  }
  expect_identical(CmxSim, mxSim)

  # Test imputing value for all NAs
  xNA <- rep(NA, length(x))
  mNA <- nimbleModel(nc, data = list(x = xNA),
         inits = list(probOcc = probOcc,
                      probDetect = probDetect))
  mNAConf <- configureMCMC(mNA)
  mNAConf$addMonitors('x')
  mNA_MCMC <- buildMCMC(mNAConf)
  cmNA <- compileNimble(mNA, mNA_MCMC)
  set.seed(0)
  cmNA$mNA_MCMC$run(10)
# Did the imputed values come back?
  expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x[1]"])))
})


test_that("Checking errors", {
### Uncompiled errors
# dOcc_ss error checks
  expect_error(
    dOcc_s(x = c(0,1,0,0), probOcc = 0.4, probDetect = 0.5, len = 3)
  )

# dOcc_sv error checks
  expect_error(
    dOcc_s(x = c(0,1,0,0), probOcc = 0.1, probDetect = c(0.9, 0.9, 0.4), len = 5)
  )
  expect_error(
    dOcc_v(x = c(0,1,0,0), probOcc = 0.1, probDetect = c(0.9, 0.9, 0.4), len = 5)
  )


  ### Compiled errors
  # we must wrap the call to avoid the error that
  # dOcc_[s,v] has no setup code but buildDerivs=TRUE so can't be compiled alone
  call_dOcc_s <- nimbleFunction(
    run = function(x = double(1), probOcc=double(),
                   probDetect=double(),
                   len = integer(0),
                   log = integer(0, default = 0)) {
      return(dOcc_s(x, probOcc, probDetect, len, log))
      returnType(double())
    }
  )
  call_dOcc_v <- nimbleFunction(
    run = function(x = double(1), probOcc=double(0),
                   probDetect=double(1),
                   len = integer(0),
                   log = integer(0, default = 0)) {
      return(dOcc_v(x, probOcc, probDetect, len, log))
      returnType(double())
    }
  )
  CdOcc_s <- compileNimble(call_dOcc_s)
  CdOcc_v <- compileNimble(call_dOcc_v)

  expect_error(
    CdOcc_s(x = c(0,1,0,0), probOcc = 0.4, probDetect = 0.5, len = 3)
  )
  expect_error(
    CdOcc_v(x = c(0,1,0,0), probOcc = 0.4, probDetect = c(0.5,0.5,0.5,0.6), len = 3)
  )

  # This should probably be set up to error:
    # expect_error(
    #   CdOcc_s(x = c(0,1,0,0), probOcc = 0.4, probDetect = c(0.5,0.5), len = 4)
    # )
  expect_error(
    CdOcc_v(x = c(0,1,0,0), probOcc = 0.4, probDetect = 0.5, len = 4)
  )

})

