# Test the Cormack-Jolly-Seber distribution nimbleFunction.

# -----------------------------------------------------------------------------
# 0. Load

# Packages
library(testthat)
library(nimble)
# Source the file
source("R/dHMM.R")
# Set the context for testthat
context("Testing dHMM-related functions.")

# -----------------------------------------------------------------------------
# 1. Test dHMM
test_that("dHMM works", {
  # Let's do a very simple example with seeds.
  # Seeds have 3 true states (s = 3):
  #     1. Alive but not sprouted
  #     2. Sprouted (visible as plant)
  #     3. Dead
  # but only 2 possible observational states (o = 2):
  #     1. Not sprouted
  #     2. Sprouted
  # (this example is going to be a bit silly because I don't think
  # capture/recapture makes sense for plants)

  # len == length(x) = t
  len <- 5
  x1 <- c(1, 1, 1, 2, 1)
  x2 <- c(1, 1, 1, 1, 1)
  # length(init) == s
  init <- c(0.4, 0.2, 0.4)

  # Z is observation probabilities, dim o x s where o is possible system states
  # (domain of x). Z says, for each true system state (row), what is the
  # corresponding probability of observing each response (col).
  Z <- t(array(
         c(1, 0.2, 1,
           0, 0.8, 0),
         c(3, 2))) # I'm doing the transpose to specify by rows, but remember that nxxx is backwards

  # Tt is transition probabilities, s x s.
  Tt <- t(array(
          c(0.6, 0.3, 0.1,
            0, 0.7, 0.3,
            0, 0, 1),
          c(3,3)))


  correctProbX1 <- 1
  correctProbX2 <- 1
  pi1 <- init
  pi2 <- init
  for (i in 1:len) {
    stateprob1 <- pi1 * Z[x1[i],]
    stateprob2 <- pi2 * Z[x2[i],]

    sumZ1 <- sum(stateprob1)
    sumZ2 <- sum(stateprob2)

    correctProbX1 <- correctProbX1 * sumZ1
    correctProbX2 <- correctProbX2 * sumZ2

    pi1 <- (Tt[,] %*% asCol(stateprob1) / sumZ1)[ ,1]
    pi2 <- (Tt[,] %*% asCol(stateprob2) / sumZ2)[ ,1]
  }


  probX1 <- dHMM(x = x1, init = init,
                 Z = Z, T = Tt,
                 len = len, log = F)
  probX2 <- dHMM(x = x2, init = init,
                 Z = Z, T = Tt,
                 len = len, log = F)

  expect_equal(probX1, correctProbX1)
  expect_equal(probX2, correctProbX2)

  lProbX1 <- dHMM(x = x1, init = init,
                  Z = Z, T = Tt,
                  len = len, log = T)
  lProbX2 <- dHMM(x = x2, init = init,
                  Z = Z, T = Tt,
                  len = len, log = T)

  expect_equal(lProbX1, log(correctProbX1))
  expect_equal(lProbX2, log(correctProbX2))


  CdHMM <- compileNimble(dHMM)
  CprobX1 <- CdHMM(x = x1, init = init,
                   Z = Z, T = Tt,
                   len = len, log = F)
  expect_equal(CprobX1, probX1)

  ClProbX1 <- CdHMM(x = x1, init = init,
                    Z = Z, T = Tt,
                    len = len, log = T)
  expect_equal(ClProbX1, lProbX1)


  ## r function works
  set.seed(1234)
  oneSim <- rHMM(1, init, Z, Tt, len = 5)
  expect_equal(oneSim, c(2, 2, 1, 2, 2))

  ## Use in model
  nc <- nimbleCode({
    # x[1:4] ~ dCJSvv(probSurvive[1:4], probCapture[1:4], len = 4)
    # for (i in 1:4) {e
    #   probSurvive[i] ~ dunif(0,1)
    #   probCapture[i] ~ dunif(0,1)
    # }
  })
  m <- nimbleModel(nc, data = list(x = x),
                   inits = list(probSurvive = probSurvive,
                                probCapture = probCapture))

})


test_that("dHMMo works", {
  # Let's do a very simple example with seeds.
  # Seeds have 3 true states (s = 3):
  #     1. Alive but not sprouted
  #     2. Sprouted (visible as plant)
  #     3. Dead
  # but only 2 possible observational states (o = 2):
  #     1. Not sprouted
  #     2. Sprouted
  # (this example is going to be a bit silly because I don't think
  # capture/recapture makes sense for plants)

  # len == length(x) = t
  len <- 5
  x1 <- c(1, 1, 1, 2, 1)
  x2 <- c(1, 1, 1, 1, 1)
  # length(init) == s
  init <- c(0.4, 0.2, 0.4)

  # Z is observation probabilities, dim o x s where o is possible system states
  # (domain of x). Z says, for each true system state (row), what is the
  # corresponding probability of observing each response (col).
  Z <- array(
         c(1, 0, 0.2, 0.8, 1, 0,
           0.9, 0.1, 0.2, 0.8, 1, 0,
           1, 0, 0.2, 0.8, 1, 0,
           1, 0, 0.2, 0.8, 1, 0,
           0.95, 0.05, 0.2, 0.8, 0.5, 0.5),
         c(2, 3, 5))

  # Tt is transition probabilities, s x s.
  Tt <- t(array(
          c(0.6, 0.3, 0.1,
            0, 0.7, 0.3,
            0, 0, 1),
          c(3,3)))


  correctProbX1 <- 1
  correctProbX2 <- 1
  pi1 <- init
  pi2 <- init
  for (i in 1:len) {
    stateprob1 <- pi1 * Z[x1[i],,i]
    stateprob2 <- pi2 * Z[x2[i],,i]

    sumZ1 <- sum(stateprob1)
    sumZ2 <- sum(stateprob2)

    correctProbX1 <- correctProbX1 * sumZ1
    correctProbX2 <- correctProbX2 * sumZ2

    pi1 <- (Tt[,] %*% asCol(stateprob1) / sumZ1)[ ,1]
    pi2 <- (Tt[,] %*% asCol(stateprob2) / sumZ2)[ ,1]
  }


  probX1 <- dHMMo(x = x1, init = init,
                  Z = Z, T = Tt,
                  len = len, log = F)
  probX2 <- dHMMo(x = x2, init = init,
                  Z = Z, T = Tt,
                  len = len, log = F)

  expect_equal(probX1, correctProbX1)
  expect_equal(probX2, correctProbX2)

  lProbX1 <- dHMMo(x = x1, init = init,
                   Z = Z, T = Tt,
                   len = len, log = T)
  lProbX2 <- dHMMo(x = x2, init = init,
                   Z = Z, T = Tt,
                   len = len, log = T)

  expect_equal(lProbX1, log(correctProbX1))
  expect_equal(lProbX2, log(correctProbX2))


  CdHMMo <- compileNimble(dHMMo)
  CprobX1 <- CdHMMo(x = x1, init = init,
                    Z = Z, T = Tt,
                    len = len, log = F)
  expect_equal(CprobX1, probX1)

  ClProbX1 <- CdHMMo(x = x1, init = init,
                     Z = Z, T = Tt,
                     len = len, log = T)
  expect_equal(ClProbX1, lProbX1)

})




test_that("dHMM errors where expected", {

# Start with good stuff and break it one by one
  len <- 5
  x <- c(1, 1, 1, 2, 1)
  init <- c(0.4, 0.2, 0.4)

  Z <- t(matrix(
         c(1, 0.2, 1,
           0, 0.8, 0),
         nrow = length(init)))

  Tt <- t(matrix(
          c(0.6, 0.3, 0.1,
            0, 0.7, 0.3,
            0, 0, 1),
          ncol = length(init)))
  badT <- t(matrix(
          c(0.6, 0.3, 0.1,
            0, 0.7, 0.3),
          ncol = length(init)))
  badZ <- t(matrix(
          c(0.6, 0.3, 0.1,
            0, 0.7, 0.3),
          ncol = length(init)))





  probX <- expect_error(
              dHMM(x = x, init = init,
                   Z = Z, T = Tt,
                   len = 4, log = F))
  probX <- expect_error(
              dHMM(x = x, init = init,
                   Z = Z, T = badT,
                   len = len, log = F))
  probX <- expect_error(
              dHMM(x = x, init = init,
                   Z = badZ, T = Tt,
                   len = len, log = F))




})



