# Test the Cormack-Jolly-Seber distribution nimbleFunction.

# -----------------------------------------------------------------------------
# 0. Load

# Packages
library(testthat)
library(nimble)
# Source the file
source("R/dCJS.R")
# Set the context for testthat
context("Testing dCJS-related functions.")

# -----------------------------------------------------------------------------
# 1. Test dCJSss
#    dCJSss is used in the case that probSurvive and probCapture are both
#    scalar values.
test_that("dCJSss works",
          {
      # Uncompiled calculation
          x <- c(0, 1, 0, 0)
          probSurvive <- 0.6
          probCapture <- 0.4
          probX <- dCJSss(x, probSurvive, probCapture)
      # Manually calculate the correct answer
          correctProbX <- probSurvive * (1 - probCapture) *
            probSurvive * (probCapture) *
            (probSurvive^2 * (1 - probCapture)^2 +
               probSurvive * (1 - probCapture) * (1 - probSurvive) +
               (1 - probSurvive))

          expect_equal(probX, correctProbX)

      # Uncompiled log probability
          lProbX <- dCJSss(x, probSurvive, probCapture, log = TRUE)
          lCorrectProbX <- log(correctProbX)
          expect_equal(lProbX, lCorrectProbX)

      # Compilation and compiled calculations
          CdCJSss <- compileNimble(dCJSss)
          CprobX <- CdCJSss(x, probSurvive, probCapture)
          expect_equal(CprobX, probX)

          ClProbX <- CdCJSss(x, probSurvive, probCapture, log = TRUE)
          expect_equal(ClProbX, lProbX)

      # random function works
          set.seed(2468)
          oneSim <- rCJSss(1, probSurvive, probCapture, len = 4)
          expect_equal(oneSim, c(0, 1, 0, 0))

      # Use in Nimble model
          nc <- nimbleCode({
            x[1:4] ~ dCJSss(probSurvive, probCapture, len = 4)
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

      # Simulate some data to test random generation
          set.seed(2468)
          cm$simulate('x')
          expect_equal(cm$x, c(0, 1, 0, 0))

      # Test imputing value for all NAs
# TODO: how do I do this for only a few NAs? The problem is right now all of x
# is a single node (which might be necessary since they're not independent?)
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
        })

# -----------------------------------------------------------------------------
# 2. Test dCJSsv
#    dCJSsv is used in the case where survival probability is a scalar and
#    capture probability is a vector.
test_that("dCJSsv works",
          {
      # Uncompiled calculation
          x <- c(0, 1, 0, 0)
          probSurvive <- 0.6
          probCapture <- c(0.25, 0.6, 0.4, 0.8)
          probX <- dCJSsv(x, probSurvive, probCapture)
      # Manually calculate the correct answer
          correctProbX <- probSurvive * (1 - probCapture[1]) *
            probSurvive * (probCapture[2]) *
            (probSurvive^2 * (1 - probCapture[3]) * (1 - probCapture[4]) +
               probSurvive * (1 - probCapture[3]) * (1 - probSurvive) +
               (1 - probSurvive))

          expect_equal(probX, correctProbX)

      # Uncompiled log probability
          lProbX <- dCJSsv(x, probSurvive, probCapture, log = TRUE)
          lCorrectProbX <- log(correctProbX)
          expect_equal(lProbX, lCorrectProbX)

      # Compilation and compiled calculations
          CdCJSsv <- compileNimble(dCJSsv)
          CprobX <- CdCJSsv(x, probSurvive, probCapture)
          expect_equal(CprobX, probX)

          ClProbX <- CdCJSsv(x, probSurvive, probCapture, log = TRUE)
          expect_equal(ClProbX, lProbX)

      # random function works
          set.seed(2)
          oneSim <- rCJSsv(1, probSurvive, probCapture, len = 4)
          expect_equal(oneSim, c(0, 1, 0, 0))

      # Use in Nimble model
          nc <- nimbleCode({
            x[1:4] ~ dCJSsv(probSurvive, probCapture[1:4], len = 4)
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
        })

# -----------------------------------------------------------------------------
# 3. Test dCJSvs
#    dCJSvs is used in the case where survival probability is a vector and
#    capture probability is a scalar.
test_that("dCJSvs works",
          {
      # Uncompiled calculation
          x <- c(0, 1, 0, 0)
          probSurvive <- c(0.8, 0.45, 0.4, 0.7)
          probCapture <- 0.6
          probX <- dCJSvs(x, probSurvive, probCapture)
      # Manually calculate the correct answer
          correctProbX <- probSurvive[1] * (1 - probCapture) *
            probSurvive[2] * probCapture *
            (probSurvive[3] * probSurvive[4] * (1 - probCapture)^2 +
               probSurvive[3] * (1 - probCapture) * (1 - probSurvive[4]) +
               (1 - probSurvive[3]))

          expect_equal(probX, correctProbX)

      # Uncompiled log probability
          lProbX <- dCJSvs(x, probSurvive, probCapture, log = TRUE)
          lCorrectProbX <- log(correctProbX)
          expect_equal(lProbX, lCorrectProbX)

      # Compilation and compiled calculations
          CdCJSvs <- compileNimble(dCJSvs)
          CprobX <- CdCJSvs(x, probSurvive, probCapture)
          expect_equal(CprobX, probX)

          ClProbX <- CdCJSvs(x, probSurvive, probCapture, log = TRUE)
          expect_equal(ClProbX, lProbX)

      # random function works
          set.seed(2)
          oneSim <- rCJSvs(1, probSurvive, probCapture, len = 4)
          expect_equal(oneSim, c(0, 1, 0, 0))

      # Use in Nimble model
          nc <- nimbleCode({
            x[1:4] ~ dCJSvs(probSurvive[1:4], probCapture, len = 4)
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
# TODO: how do I do this for only a few NAs? The problem is right now all of x
# is a single node (which might be necessary since they're not independent?)
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
        })

# -----------------------------------------------------------------------------
# 4. Test dCJSvv

test_that("dCJSvv works",
      {
      ## Uncompiled calculation
          x <- c(0,1,0,0)
          probSurvive <- c(0.6, 0.5, 0.4, 0.55)
          probCapture <- c(0.45, 0.5, 0.55, 0.6)
          len <- 4
          probX <- dCJSvv(x, probSurvive, probCapture, len)

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
          lProbX <- dCJSvv(x, probSurvive, probCapture, log = TRUE)
          lCorrectProbX <- log(correctProbX)
          expect_equal(lProbX, lCorrectProbX)

          ## Compiles
          CdCJSvv <- compileNimble(dCJSvv)
          CprobX <- CdCJSvv(x, probSurvive, probCapture)
          expect_equal(CprobX, probX)

          ClProbX <- CdCJSvv(x, probSurvive, probCapture, len = 4, log = TRUE)
          expect_equal(ClProbX, lProbX)

          ## r function works
          set.seed(2468)
          oneSim <- rCJSvv(1, probSurvive, probCapture, len = 4)
          expect_equal(oneSim, c(1, 0, 0, 0))

          ## Use in model
          nc <- nimbleCode({
            x[1:4] ~ dCJSvv(probSurvive[1:4], probCapture[1:4], len = 4)
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
# TODO: how do I do this for only a few NAs? The problem is right now all of x
# is a single node (which might be necessary since they're not independent?)
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
  })

test_that("dCJS errors", {

### Uncompiled errors
# dCJSss error checks
  expect_error(
    dCJSss(x = c(0,1,0,0), probCapture = 0.4, probSurvive = 0.5, len = 3)
  )

# dCJSsv error checks
  expect_error(
    dCJSsv(x = c(0,1,0,0), probCapture = 0.1, probSurvive = c(0.9, 0.9, 0.4, 0.4), len = 5)
  )
  expect_error(
    dCJSsv(x = c(0,1,0,0), probCapture = 0.1, probSurvive = c(0.9, 0.9, 0.4))
  )

# dCJSvs error checks
  expect_error(
    dCJSsv(x = c(0,1,0,0), probCapture = c(0.9, 0.9, 0.4, 0.4), probSurvive = 0.1, len = 5)
  )
  expect_error(
    dCJSsv(x = c(0,1,0,0), probCapture = c(0.9, 0.9, 0.4), probSurvive = 0.8)
  )

# dCJSvv error checks
  expect_error(
    dCJSvv(x = c(0,1,0,0), probCapture = c(0,1,0.3,0.3), probSurvive = c(0.9, 0.9))
  )
  expect_error(
    dCJSvv(x = c(0,1,0,0), probCapture = c(0,1), probSurvive = c(0.9, 0.9, 0.1, 0.1))
  )
  expect_error(
    dCJSvv(x = c(0,1,0,0), probCapture = c(0,1,0,0),
            probSurvive = c(0.9, 0.9, 0.1, 0.1), len = 2)
  )

### Compiled errors
  CdCJSss <- compileNimble(dCJSss)
  CdCJSsv <- compileNimble(dCJSsv)
  CdCJSvs <- compileNimble(dCJSvs)
  CdCJSvv <- compileNimble(dCJSvv)

  expect_error(
    CdCJSss(x = c(0,1,0,0), probCapture = 0.4, probSurvive = 0.5, len = 3)
  )

# dCJSsv error checks
  expect_error(
    CdCJSsv(x = c(0,1,0,0), probCapture = 0.1, probSurvive = c(0.9, 0.9, 0.4, 0.4), len = 5)
  )
  expect_error(
    CdCJSsv(x = c(0,1,0,0), probCapture = 0.1, probSurvive = c(0.9, 0.9, 0.4))
  )

# dCJSvs error checks
  expect_error(
    CdCJSsv(x = c(0,1,0,0), probCapture = c(0.9, 0.9, 0.4, 0.4), probSurvive = 0.1, len = 5)
  )
  expect_error(
    CdCJSsv(x = c(0,1,0,0), probCapture = c(0.9, 0.9, 0.4), probSurvive = 0.8)
  )

# dCJSvv error checks
  expect_error(
    CdCJSvv(x = c(0,1,0,0), probCapture = c(0,1,0.3,0.3), probSurvive = c(0.9, 0.9))
  )
  expect_error(
    CdCJSvv(x = c(0,1,0,0), probCapture = c(0,1), probSurvive = c(0.9, 0.9, 0.1, 0.1))
  )
  expect_error(
    CdCJSvv(x = c(0,1,0,0), probCapture = c(0,1,0,0),
            probSurvive = c(0.9, 0.9, 0.1, 0.1), len = 2)
  )
})

