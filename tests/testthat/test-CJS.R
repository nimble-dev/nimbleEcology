# Test the Cormack-Jolly-Seber distribution nimbleFunction.

# -----------------------------------------------------------------------------
# 0. Load
# Set the context for testthat
context("Testing dCJS-related functions.")

# -----------------------------------------------------------------------------
# 1. Test dCJS_ss
#    dCJS_ss is used in the case that probSurvive and probCapture are both
#    scalar values.
test_that("dCJS_ss works",
          {
      # Uncompiled calculation
          x <- c(0, 1, 0, 0)
          probSurvive <- 0.6
          probCapture <- 0.4
          probX <- dCJS_ss(x, probSurvive, probCapture)
      # Manually calculate the correct answer
          correctProbX <- probSurvive * (1 - probCapture) *
            probSurvive * (probCapture) *
            (probSurvive^2 * (1 - probCapture)^2 +
               probSurvive * (1 - probCapture) * (1 - probSurvive) +
               (1 - probSurvive))

          expect_equal(probX, correctProbX)

      # Uncompiled log probability
          lProbX <- dCJS_ss(x, probSurvive, probCapture, log = TRUE)
          lCorrectProbX <- log(correctProbX)
          expect_equal(lProbX, lCorrectProbX)

      # Compilation and compiled calculations
          CdCJS_ss <- compileNimble(dCJS_ss)
          CprobX <- CdCJS_ss(x, probSurvive, probCapture)
          expect_equal(CprobX, probX)

          ClProbX <- CdCJS_ss(x, probSurvive, probCapture, log = TRUE)
          expect_equal(ClProbX, lProbX)

      # random function works
          set.seed(2468)
          oneSim <- rCJS_ss(1, probSurvive, probCapture, len = 4)
          expect_equal(oneSim, c(0, 1, 0, 0))

      # Use in Nimble model
          nc <- nimbleCode({
            x[1:4] ~ dCJS_ss(probSurvive, probCapture, len = 4)
            probSurvive ~ dunif(0,1)
            probCapture ~ dunif(0,1)
          })
          m <- nimbleModel(nc, data = list(x = x),
                           inits = list(probSurvive = probSurvive,
                                        probCapture = probCapture))
          m$calculate()
          MlProbX <- m$getLogProb("x")
          expect_equal(MlProbX, lProbX)

      # Compiled model
          cm <- compileNimble(m)
          cm$calculate()
          CMlProbX <- cm$getLogProb("x")
          expect_equal(CMlProbX, lProbX)

      # Test imputing value for all NAs
          xNA <- c(NA, NA, NA, NA)
          mNA <- nimbleModel(nc, data = list(x = xNA),
                 inits = list(probSurvive = probSurvive,
                              probCapture = probCapture))
          mNAConf <- configureMCMC(mNA)
          mNAConf$addMonitors('x')
          mNA_MCMC <- buildMCMC(mNAConf)
          cmNA <- compileNimble(mNA, mNA_MCMC)

          set.seed(0)
          cmNA$mNA_MCMC$run(10)

      # Did the imputed values come back?
          expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x[1]"])))

      # Test simulation code
          set.seed(1)
          nSim <- 10
          xSim <- array(NA, dim = c(nSim, length(x)))
          for(i in 1:nSim)
            xSim[i,] <- rCJS_ss(1, probSurvive, probCapture, len = length(x))
          set.seed(1)
          CrCJS_ss <- compileNimble(rCJS_ss)
          CxSim <- array(NA, dim = c(nSim, length(x)))
          for(i in 1:nSim)
            CxSim[i,] <- CrCJS_ss(1, probSurvive, probCapture, len = length(x))
          expect_identical(xSim, CxSim)

          simNodes <- m$getDependencies(c('probSurvive', 'probCapture'), self = FALSE)
          mxSim <- array(NA, dim = c(nSim, length(x)))
          set.seed(1)
          for(i in 1:nSim) {
            m$simulate(simNodes, includeData = TRUE)
            mxSim[i,] <- m$x
          }
          expect_identical(mxSim, xSim)

          CmxSim <- array(NA, dim = c(nSim, length(x)))
          set.seed(1)
          for(i in 1:nSim) {
            cm$simulate(simNodes, includeData = TRUE)
            CmxSim[i,] <- cm$x
          }
          expect_identical(CmxSim, mxSim)
        })

# -----------------------------------------------------------------------------
# 2. Test dCJS_sv
#    dCJS_sv is used in the case where survival probability is a scalar and
#    capture probability is a vector.
test_that("dCJS_sv works",
          {
      # Uncompiled calculation
          x <- c(0, 1, 0, 0)
          probSurvive <- 0.6
          probCapture <- c(0.25, 0.6, 0.4, 0.8)
          probX <- dCJS_sv(x, probSurvive, probCapture)
      # Manually calculate the correct answer
          correctProbX <- probSurvive * (1 - probCapture[1]) *
            probSurvive * (probCapture[2]) *
            (probSurvive^2 * (1 - probCapture[3]) * (1 - probCapture[4]) +
               probSurvive * (1 - probCapture[3]) * (1 - probSurvive) +
               (1 - probSurvive))

          expect_equal(probX, correctProbX)

      # Uncompiled log probability
          lProbX <- dCJS_sv(x, probSurvive, probCapture, log = TRUE)
          lCorrectProbX <- log(correctProbX)
          expect_equal(lProbX, lCorrectProbX)

      # Compilation and compiled calculations
          CdCJS_sv <- compileNimble(dCJS_sv)
          CprobX <- CdCJS_sv(x, probSurvive, probCapture)
          expect_equal(CprobX, probX)

          ClProbX <- CdCJS_sv(x, probSurvive, probCapture, log = TRUE)
          expect_equal(ClProbX, lProbX)

      # random function works
          set.seed(2)
          oneSim <- rCJS_sv(1, probSurvive, probCapture, len = 4)
          expect_equal(oneSim, c(0, 1, 0, 0))

      # Use in Nimble model
          nc <- nimbleCode({
            x[1:4] ~ dCJS_sv(probSurvive, probCapture[1:4], len = 4)
            probSurvive ~ dunif(0,1)
            for (i in 1:4) {
              probCapture[i] ~ dunif(0,1)
            }
          })
          m <- nimbleModel(nc, data = list(x = x),
                           inits = list(probSurvive = probSurvive,
                                        probCapture = probCapture))
          m$calculate()
          MlProbX <- m$getLogProb("x")
          expect_equal(MlProbX, lProbX)

      # Compiled model
          cm <- compileNimble(m)
          cm$calculate()
          CMlProbX <- cm$getLogProb("x")
          expect_equal(CMlProbX, lProbX)

      # Simulate some data to test random generation
          set.seed(2468)
          cm$simulate('x')
          expect_equal(cm$x, c(0, 1, 0, 0))

      # Test imputing value for all NAs
          xNA <- c(NA, NA, NA, NA)
          mNA <- nimbleModel(nc, data = list(x = xNA),
                 inits = list(probSurvive = probSurvive,
                              probCapture = probCapture))
          mNAConf <- configureMCMC(mNA)
          mNAConf$addMonitors('x')
          mNA_MCMC <- buildMCMC(mNAConf)
          cmNA <- compileNimble(mNA, mNA_MCMC)

          set.seed(10)
          cmNA$mNA_MCMC$run(5)

      # Did the imputed values come back?
          expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x[1]"])))

      # Test simulation code
          set.seed(1)
          nSim <- 10
          xSim <- array(NA, dim = c(nSim, length(x)))
          for(i in 1:nSim)
            xSim[i,] <- rCJS_sv(1, probSurvive, probCapture, len = length(x))
          set.seed(1)
          CrCJS_sv <- compileNimble(rCJS_sv)
          CxSim <- array(NA, dim = c(nSim, length(x)))
          for(i in 1:nSim)
            CxSim[i,] <- CrCJS_sv(1, probSurvive, probCapture, len = length(x))
          expect_identical(xSim, CxSim)

          simNodes <- m$getDependencies(c('probSurvive', 'probCapture'), self = FALSE)
          mxSim <- array(NA, dim = c(nSim, length(x)))
          set.seed(1)
          for(i in 1:nSim) {
            m$simulate(simNodes, includeData = TRUE)
            mxSim[i,] <- m$x
          }
          expect_identical(mxSim, xSim)

          CmxSim <- array(NA, dim = c(nSim, length(x)))
          set.seed(1)
          for(i in 1:nSim) {
            cm$simulate(simNodes, includeData = TRUE)
            CmxSim[i,] <- cm$x
          }
          expect_identical(CmxSim, mxSim)
        })

# -----------------------------------------------------------------------------
# 3. Test dCJS_vs
#    dCJS_vs is used in the case where survival probability is a vector and
#    capture probability is a scalar.
test_that("dCJS_vs works",
          {
      # Uncompiled calculation
          x <- c(0, 1, 0, 0)
          probSurvive <- c(0.8, 0.45, 0.4, 0.7)
          probCapture <- 0.6
          probX <- dCJS_vs(x, probSurvive, probCapture)
      # Manually calculate the correct answer
          correctProbX <- probSurvive[1] * (1 - probCapture) *
            probSurvive[2] * probCapture *
            (probSurvive[3] * probSurvive[4] * (1 - probCapture)^2 +
               probSurvive[3] * (1 - probCapture) * (1 - probSurvive[4]) +
               (1 - probSurvive[3]))

          expect_equal(probX, correctProbX)

      # Uncompiled log probability
          lProbX <- dCJS_vs(x, probSurvive, probCapture, log = TRUE)
          lCorrectProbX <- log(correctProbX)
          expect_equal(lProbX, lCorrectProbX)

      # Compilation and compiled calculations
          CdCJS_vs <- compileNimble(dCJS_vs)
          CprobX <- CdCJS_vs(x, probSurvive, probCapture)
          expect_equal(CprobX, probX)

          ClProbX <- CdCJS_vs(x, probSurvive, probCapture, log = TRUE)
          expect_equal(ClProbX, lProbX)

      # random function works
          set.seed(2)
          oneSim <- rCJS_vs(1, probSurvive, probCapture, len = 4)
          expect_equal(oneSim, c(0, 1, 0, 0))

      # Use in Nimble model
          nc <- nimbleCode({
            x[1:4] ~ dCJS_vs(probSurvive[1:4], probCapture, len = 4)
            probCapture ~ dunif(0,1)
            for (i in 1:4) {
              probSurvive[i] ~ dunif(0,1)
            }
          })
          m <- nimbleModel(nc, data = list(x = x),
                           inits = list(probSurvive = probSurvive,
                                        probCapture = probCapture))
          m$calculate()
          MlProbX <- m$getLogProb("x")
          expect_equal(MlProbX, lProbX)

      # Compiled model
          cm <- compileNimble(m)
          cm$calculate()
          CMlProbX <- cm$getLogProb("x")
          expect_equal(CMlProbX, lProbX)

      # Simulate some data to test random generation
          set.seed(2468)
          cm$simulate('x')
          expect_equal(cm$x, c(0, 1, 0, 0))

      # Test imputing value for all NAs
          xNA <- c(NA, NA, NA, NA)
          mNA <- nimbleModel(nc, data = list(x = xNA),
                 inits = list(probSurvive = probSurvive,
                              probCapture = probCapture))
          mNAConf <- configureMCMC(mNA)
          mNAConf$addMonitors('x')
          mNA_MCMC <- buildMCMC(mNAConf)
          cmNA <- compileNimble(mNA, mNA_MCMC)

          set.seed(5)
          cmNA$mNA_MCMC$run(10)

      # Did the imputed values come back?
          expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x[1]"])))

      # Test simulation code
          set.seed(1)
          nSim <- 10
          xSim <- array(NA, dim = c(nSim, length(x)))
          for(i in 1:nSim)
            xSim[i,] <- rCJS_vs(1, probSurvive, probCapture, len = length(x))
          set.seed(1)
          CrCJS_vs <- compileNimble(rCJS_vs)
          CxSim <- array(NA, dim = c(nSim, length(x)))
          for(i in 1:nSim)
            CxSim[i,] <- CrCJS_vs(1, probSurvive, probCapture, len = length(x))
          expect_identical(xSim, CxSim)

          simNodes <- m$getDependencies(c('probSurvive', 'probCapture'), self = FALSE)
          mxSim <- array(NA, dim = c(nSim, length(x)))
          set.seed(1)
          for(i in 1:nSim) {
            m$simulate(simNodes, includeData = TRUE)
            mxSim[i,] <- m$x
          }
          expect_identical(mxSim, xSim)

          CmxSim <- array(NA, dim = c(nSim, length(x)))
          set.seed(1)
          for(i in 1:nSim) {
            cm$simulate(simNodes, includeData = TRUE)
            CmxSim[i,] <- cm$x
          }
          expect_identical(CmxSim, mxSim)
        })

# -----------------------------------------------------------------------------
# 4. Test dCJS_vv

test_that("dCJS_vv works",
      {
      ## Uncompiled calculation
          x <- c(0,1,0,0)
          probSurvive <- c(0.6, 0.5, 0.4, 0.55)
          probCapture <- c(0.45, 0.5, 0.55, 0.6)
          len <- 4
          probX <- dCJS_vv(x, probSurvive, probCapture, len)

          correctProbX <-
            probSurvive[1] * (1 - probCapture[1]) *
            probSurvive[2] * (probCapture[2]) *
            ((probSurvive[3] * (1 - probCapture[3]) *
              probSurvive[4] * (1 - probCapture[4])) +
             (probSurvive[3] * (1 - probCapture[3]) *
              (1 - probSurvive[4])) +
             (1 - probSurvive[3]))

          expect_equal(probX, correctProbX)

          ## log Prob
          lProbX <- dCJS_vv(x, probSurvive, probCapture, log = TRUE)
          lCorrectProbX <- log(correctProbX)
          expect_equal(lProbX, lCorrectProbX)

          ## Compiles
          CdCJS_vv <- compileNimble(dCJS_vv)
          CprobX <- CdCJS_vv(x, probSurvive, probCapture)
          expect_equal(CprobX, probX)

          ClProbX <- CdCJS_vv(x, probSurvive, probCapture, len = 4, log = TRUE)
          expect_equal(ClProbX, lProbX)

          ## r function works
          set.seed(2468)
          oneSim <- rCJS_vv(1, probSurvive, probCapture, len = 4)
          expect_equal(oneSim, c(1, 0, 0, 0))

          ## Use in model
          nc <- nimbleCode({
            x[1:4] ~ dCJS_vv(probSurvive[1:4], probCapture[1:4], len = 4)
            for (i in 1:4) {
              probSurvive[i] ~ dunif(0,1)
              probCapture[i] ~ dunif(0,1)
            }
          })
          m <- nimbleModel(nc, data = list(x = x),
                           inits = list(probSurvive = probSurvive,
                                        probCapture = probCapture))
          m$calculate()
          MlProbX <- m$getLogProb("x")
          expect_equal(MlProbX, lProbX)

          cm <- compileNimble(m)
          cm$calculate()
          CMlProbX <- cm$getLogProb("x")
          expect_equal(CMlProbX, lProbX)

          set.seed(2468)
          cm$simulate('x')
          expect_equal(cm$x, c(0, 1, 0, 0))

      # Test imputing value for all NAs
          xNA <- c(NA, NA, NA, NA)
          mNA <- nimbleModel(nc, data = list(x = xNA),
                 inits = list(probSurvive = probSurvive,
                              probCapture = probCapture))
          mNAConf <- configureMCMC(mNA)
          mNAConf$addMonitors('x')
          mNA_MCMC <- buildMCMC(mNAConf)
          cmNA <- compileNimble(mNA, mNA_MCMC)

          set.seed(5)
          cmNA$mNA_MCMC$run(10)

      # Did the imputed values come back?
          expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x[1]"])))

      # Test simulation code
          set.seed(1)
          nSim <- 10
          xSim <- array(NA, dim = c(nSim, length(x)))
          for(i in 1:nSim)
            xSim[i,] <- rCJS_vv(1, probSurvive, probCapture, len = length(x))
          set.seed(1)
          CrCJS_vv <- compileNimble(rCJS_vv)
          CxSim <- array(NA, dim = c(nSim, length(x)))
          for(i in 1:nSim)
            CxSim[i,] <- CrCJS_vv(1, probSurvive, probCapture, len = length(x))
          expect_identical(xSim, CxSim)

          simNodes <- m$getDependencies(c('probSurvive', 'probCapture'), self = FALSE)
          mxSim <- array(NA, dim = c(nSim, length(x)))
          set.seed(1)
          for(i in 1:nSim) {
            m$simulate(simNodes, includeData = TRUE)
            mxSim[i,] <- m$x
          }
          expect_identical(mxSim, xSim)

          CmxSim <- array(NA, dim = c(nSim, length(x)))
          set.seed(1)
          for(i in 1:nSim) {
            cm$simulate(simNodes, includeData = TRUE)
            CmxSim[i,] <- cm$x
          }
          expect_identical(CmxSim, mxSim)
  })

test_that("dCJS errors", {

### Uncompiled errors
# dCJS_ss error checks
  expect_error(
    dCJS_ss(x = c(0,1,0,0), probCapture = 0.4, probSurvive = 0.5, len = 3)
  )

# dCJS_sv error checks
  expect_error(
    dCJS_sv(x = c(0,1,0,0), probCapture = 0.1, probSurvive = c(0.9, 0.9, 0.4, 0.4), len = 5)
  )
  expect_error(
    dCJS_sv(x = c(0,1,0,0), probCapture = 0.1, probSurvive = c(0.9, 0.9, 0.4))
  )

# dCJS_vs error checks
  expect_error(
    dCJS_sv(x = c(0,1,0,0), probCapture = c(0.9, 0.9, 0.4, 0.4), probSurvive = 0.1, len = 5)
  )
  expect_error(
    dCJS_sv(x = c(0,1,0,0), probCapture = c(0.9, 0.9, 0.4), probSurvive = 0.8)
  )

# dCJS_vv error checks
  expect_error(
    dCJS_vv(x = c(0,1,0,0), probCapture = c(0,1,0.3,0.3), probSurvive = c(0.9, 0.9))
  )
  expect_error(
    dCJS_vv(x = c(0,1,0,0), probCapture = c(0,1), probSurvive = c(0.9, 0.9, 0.1, 0.1))
  )
  expect_error(
    dCJS_vv(x = c(0,1,0,0), probCapture = c(0,1,0,0),
            probSurvive = c(0.9, 0.9, 0.1, 0.1), len = 2)
  )

### Compiled errors
  CdCJS_ss <- compileNimble(dCJS_ss)
  CdCJS_sv <- compileNimble(dCJS_sv)
  CdCJS_vs <- compileNimble(dCJS_vs)
  CdCJS_vv <- compileNimble(dCJS_vv)

  expect_error(
    CdCJS_ss(x = c(0,1,0,0), probCapture = 0.4, probSurvive = 0.5, len = 3)
  )

# dCJS_sv error checks
  expect_error(
    CdCJS_sv(x = c(0,1,0,0), probCapture = 0.1, probSurvive = c(0.9, 0.9, 0.4, 0.4), len = 5)
  )
  expect_error(
    CdCJS_sv(x = c(0,1,0,0), probCapture = 0.1, probSurvive = c(0.9, 0.9, 0.4))
  )

# dCJS_vs error checks
  expect_error(
    CdCJS_sv(x = c(0,1,0,0), probCapture = c(0.9, 0.9, 0.4, 0.4), probSurvive = 0.1, len = 5)
  )
  expect_error(
    CdCJS_sv(x = c(0,1,0,0), probCapture = c(0.9, 0.9, 0.4), probSurvive = 0.8)
  )

# dCJS_vv error checks
  expect_error(
    CdCJS_vv(x = c(0,1,0,0), probCapture = c(0,1,0.3,0.3), probSurvive = c(0.9, 0.9))
  )
  expect_error(
    CdCJS_vv(x = c(0,1,0,0), probCapture = c(0,1), probSurvive = c(0.9, 0.9, 0.1, 0.1))
  )
  expect_error(
    CdCJS_vv(x = c(0,1,0,0), probCapture = c(0,1,0,0),
            probSurvive = c(0.9, 0.9, 0.1, 0.1), len = 2)
  )
})

