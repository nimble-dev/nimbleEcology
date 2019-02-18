# tests for occupancy distribution functions

library(testthat)
library(nimble)
source("R/dOcc.R")
context("Testing dOcc-related functions.")

test_that("dOcc_s works",
          {
    x <- c(1,0,1,1,0)
    probOcc <- 0.4
    probDetect <- 0.7

    probX <- dOcc_s(x, probOcc, probDetect)
    correctProbX <-
        probOcc * probDetect^3 * (1 - probDetect)^2
    expect_equal(probX, correctProbX)

    lProbX <- dOcc_s(x, probOcc, probDetect, log = TRUE)
    lCorrectProbX <- log(correctProbX)
    expect_equal(lProbX, lCorrectProbX)

    CdOcc_s <- compileNimble(dOcc_s)
    CprobX <- CdOcc_s(x, probOcc, probDetect)
    expect_equal(CprobX, probX)

    ClProbX <- CdOcc_s(x, probOcc, probDetect, log = TRUE)
    expect_equal(ClProbX, lProbX)

    nc <- nimbleCode({
      x[1:5] ~ dOcc_s(probOcc, probDetect)
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
    expect_equal(cm$x, c(0, 1, 0, 0))

})


test_that("dOcc_v works",
          {
    x <- c(1,0,1,1,0)
    probOcc <- 0.4
    probDetect <- c(0.7, 0.3, 0.5, 0.7, 0.25)

    probX <- dOcc_v(x, probOcc, probDetect)
    correctProbX <-
        probOcc * prod(probDetect[x == 1]) *
            prod(1 - probDetect[x != 1])
    expect_equal(probX, correctProbX)

    lProbX <- dOcc_s(x, probOcc, probDetect, log = TRUE)
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




