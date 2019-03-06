# tests for occupancy distribution functions

library(testthat)
library(nimble)
source("R/dOcc.R")
context("Testing dOcc-related functions.")

# Test scalar-scalar version
test_that("dOcc_ss works",
          {
    x <- c(1,0,1,1,0)
    probOcc <- 0.4
    probDetect <- 0.7

    probX <- dOcc_ss(x, probOcc, probDetect)
    correctProbX <-
        probOcc * probDetect^3 * (1 - probDetect)^2
    expect_equal(probX, correctProbX)

    lProbX <- dOcc_ss(x, probOcc, probDetect, log = TRUE)
    lCorrectProbX <- log(correctProbX)
    expect_equal(lProbX, lCorrectProbX)

    CdOcc_ss <- compileNimble(dOcc_ss)
    CprobX <- CdOcc_ss(x, probOcc, probDetect)
    expect_equal(CprobX, probX)

    ClProbX <- CdOcc_ss(x, probOcc, probDetect, log = TRUE)
    expect_equal(ClProbX, lProbX)

    nc <- nimbleCode({
      x[1:5] ~ dOcc_ss(probOcc, probDetect)
      probDetect ~ dunif(0,1)
      probOcc ~ dunif(0,1)
    })
    m <- nimbleModel(nc, data = list(x = x),
                     inits = list(probOcc = probOcc,
                                  probDetect = probDetect))
    m$calculate()
    MlProbX <- m$getLogProb("x")
    expect_equal(MlProbX, lProbX)

    cm <- compileNimble(m, showCompilerOutput = TRUE)
    cm$calculate()
    CMlProbX <- cm$getLogProb("x")
    expect_equal(CMlProbX, lProbX)

    set.seed(2468)
    cm$simulate('x')
    expect_equal(cm$x, c(1, 0, 1, 1, 0))

})


# Test sv version
test_that("dOcc_sv works",
          {
    x <- c(1,0,1,1,0)
    probOcc <- 0.4
    probDetect <- c(0.7, 0.3, 0.5, 0.7, 0.25)

    probX <- dOcc_sv(x, probOcc, probDetect)
    correctProbX <-
        probOcc * prod(probDetect[x == 1]) *
            prod(1 - probDetect[x != 1])
    expect_equal(probX, correctProbX)

    lProbX <- dOcc_sv(x, probOcc, probDetect, log = TRUE)
    lCorrectProbX <- log(correctProbX)
    expect_equal(lProbX, lCorrectProbX)

    CdOcc_sv <- compileNimble(dOcc_sv)
    CprobX <- CdOcc_sv(x, probOcc, probDetect)
    expect_equal(CprobX, probX)

    ClProbX <- CdOcc_sv(x, probOcc, probDetect, log = TRUE)
    expect_equal(ClProbX, lProbX)

    nc <- nimbleCode({
      x[1:5] ~ dOcc_sv(probOcc, probDetect[1:5])
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

    cm <- compileNimble(m, showCompilerOutput = TRUE)
    cm$calculate()
    CMlProbX <- cm$getLogProb("x")
    expect_equal(CMlProbX, lProbX)

    set.seed(2468)
    cm$simulate('x')
    expect_equal(cm$x, c(1, 0, 1, 1, 0))

})


# Test vs version
test_that("dOcc_vs works",
          {
    x <- c(1,0,1,1,0)
    probOcc <- c(0.9, 0.4, 0.3, 0.8, 0.1)
    probDetect <- 0.7

    probX <- dOcc_vs(x, probOcc, probDetect)
    correctProbX <-
        probOcc[1] * probDetect *
        (probOcc[2] * (1 - probDetect) + (1 - probOcc[2])) *
        probOcc[3] * probDetect *
        probOcc[4] * probDetect *
        (probOcc[5] * (1 - probDetect) + (1 - probOcc[5]))


    expect_equal(probX, correctProbX)

    lProbX <- dOcc_vs(x, probOcc, probDetect, log = TRUE)
    lCorrectProbX <- log(correctProbX)
    expect_equal(lProbX, lCorrectProbX)

    CdOcc_v <- compileNimble(dOcc_v)
    CprobX <- CdOcc_v(x, probOcc, probDetect)
    expect_equal(CprobX, probX)

    ClProbX <- CdOcc_v(x, probOcc, probDetect, log = TRUE)
    expect_equal(ClProbX, lProbX)

    nc <- nimbleCode({
      x[1:5] ~ dOcc_v(probOcc, probDetect[1:5])
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

    cm <- compileNimble(m, showCompilerOutput = TRUE)
    cm$calculate()
    CMlProbX <- cm$getLogProb("x")
    expect_equal(CMlProbX, lProbX)

    set.seed(2468)
    cm$simulate('x')
    expect_equal(cm$x, c(1, 0, 1, 1, 0))

})


