# Test the Hidden Markov Model distribution nimbleFunction.

# -----------------------------------------------------------------------------
# 0. Load

# Set the context for testthat
context("Testing dHMM-related functions.")

# -----------------------------------------------------------------------------
# 1. Test dHMM, distribution for Hidden Markov Model
test_that("dHMM works", {

  # len: length of the data
  len <- 5
  # Two different data samples
  x1 <- c(1, 1, 1, 2, 1)
  x2 <- c(1, 1, 1, 1, 1)
  # length(init) == s and sum(init) == 1; initial state probabilities
  init <- c(0.4, 0.2, 0.4)

  # Z is observation probabilities, dim o x s where o is possible observation
  # states (domain of x) and s is possible true states. Z says, for each true
  # system state (row), what is the corresponding probability of observing each
  # response (col)
  Z <- t(array(
         c(1, 0.2, 1,
           0, 0.8, 0),
         c(3, 2)))

  # Tt is transition probabilities, s x s.
  Tt <- t(array(
          c(0.6, 0.3, 0.1,
            0, 0.7, 0.3,
            0, 0, 1),
          c(3,3)))

  # Calculate probabilities by hand
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

  # Calculate probabilities of x1 and x2 using dHMM
  probX1 <- dHMM(x = x1, init = init,
                 Z = Z, T = Tt,
                 len = len, log = F)
  probX2 <- dHMM(x = x2, init = init,
                 Z = Z, T = Tt,
                 len = len, log = F)

  # Compare hand-calculated and function-calculated
  expect_equal(probX1, correctProbX1)
  expect_equal(probX2, correctProbX2)

  # Repeat for log prob
  lProbX1 <- dHMM(x = x1, init = init,
                  Z = Z, T = Tt,
                  len = len, log = TRUE)
  lProbX2 <- dHMM(x = x2, init = init,
                  Z = Z, T = Tt,
                  len = len, log = TRUE)

  expect_equal(lProbX1, log(correctProbX1))
  expect_equal(lProbX2, log(correctProbX2))

  # Repeat for the compiled function
  CdHMM <- compileNimble(dHMM)
  CprobX1 <- CdHMM(x = x1, init = init,
                   Z = Z, T = Tt,
                   len = len, log = FALSE)
  expect_equal(CprobX1, probX1)

  ClProbX1 <- CdHMM(x = x1, init = init,
                    Z = Z, T = Tt,
                    len = len, log = TRUE)
  expect_equal(ClProbX1, lProbX1)

  # Create code for a nimbleModel using the distribution
  nc <- nimbleCode({
    x[1:5] ~ dHMM(init[1:3], Z = Z[1:2,1:3],
                  T = Tt[1:3, 1:3], len = 5)

    for (i in 1:3) {
      init[i] ~ dunif(0,1)

      for (j in 1:3) {
        Tt[i,j] ~ dunif(0,1)
      }

      Z[1,i] ~ dunif(0,1)
      Z[2,i] <- 1 - Z[1,i]
    }
  })

  # Create a nimbleModel using the distribution
  m <- nimbleModel(nc, data = list(x = x1),
                   inits = list(init = init, Z = Z,
                                Tt = Tt))
  # Calculate probability of x from the model
  m$calculate()
  MlProbX <- m$getLogProb("x")
  expect_equal(MlProbX, lProbX1)

  # Compile the model and re-calculate
  cm <- compileNimble(m)
  cm$calculate()
  CMlProbX <- cm$getLogProb("x")
  expect_equal(CMlProbX, lProbX1)

  # Test simulation code
  set.seed(1)
  nSim <- 10
  xSim <- array(NA, dim = c(nSim, length(x1)))
  for(i in 1:nSim)
    xSim[i,] <- rHMM(1, init, Z, Tt, len = length(x1))
  set.seed(1)
  CrHMM <- compileNimble(rHMM)
  CxSim <- array(NA, dim = c(nSim, length(x1)))
  for(i in 1:nSim)
    CxSim[i,] <- CrHMM(1, init, Z, Tt, len = length(x1))
  expect_identical(xSim, CxSim)

  simNodes <- m$getDependencies(c('init', 'Z', 'Tt'), self = FALSE)
  mxSim <- array(NA, dim = c(nSim, length(x1)))
  set.seed(1)
  for(i in 1:nSim) {
    m$simulate(simNodes, includeData = TRUE)
    mxSim[i,] <- m$x
  }
  expect_identical(mxSim, xSim)

  CmxSim <- array(NA, dim = c(nSim, length(x1)))
  set.seed(1)
  for(i in 1:nSim) {
    cm$simulate(simNodes, includeData = TRUE)
    CmxSim[i,] <- cm$x
  }
  expect_identical(CmxSim, mxSim)

  # Test imputing value for all NAs
  xNA <- re(NA, length(x))
  mNA <- nimbleModel(nc, data = list(x = xNA),
                     inits = list(init = init, Z = Z,
                                  Tt = Tt))
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
# 2. Test dHMMo, distribution for Hidden Markov Model
test_that("dHMMo works", {
  len <- 5
  x1 <- c(1, 1, 1, 2, 1)
  x2 <- c(1, 1, 1, 1, 1)
  # length(init) == s
  init <- c(0.4, 0.2, 0.4)

  # Z is observation probabilities, dim o x s x t where o is possible
  # observation states (domain of x), s is possible true system states, and t is
  # time index. Z says, for each true system state (row), what is the
  # corresponding probability of observing each response (col) at time t (3rd
  # dim)
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

  # Calculate the probabilities of x1 and x2 responses given the params
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

  # Calculate probabilities of xs using the function dHMMo
  probX1 <- dHMMo(x = x1, init = init,
                  Z = Z, T = Tt,
                  len = len, log = F)
  probX2 <- dHMMo(x = x2, init = init,
                  Z = Z, T = Tt,
                  len = len, log = F)

  # Compare function and manual
  expect_equal(probX1, correctProbX1)
  expect_equal(probX2, correctProbX2)

  # Repeat for log prob
  lProbX1 <- dHMMo(x = x1, init = init,
                   Z = Z, T = Tt,
                   len = len, log = TRUE)
  lProbX2 <- dHMMo(x = x2, init = init,
                   Z = Z, T = Tt,
                   len = len, log = TRUE)

  expect_equal(lProbX1, log(correctProbX1))
  expect_equal(lProbX2, log(correctProbX2))

  # Repeat for compiled nimbleFunction
  CdHMMo <- compileNimble(dHMMo)
  CprobX1 <- CdHMMo(x = x1, init = init,
                    Z = Z, T = Tt,
                    len = len, log = FALSE)
  expect_equal(CprobX1, probX1)

  ClProbX1 <- CdHMMo(x = x1, init = init,
                     Z = Z, T = Tt,
                     len = len, log = TRUE)
  expect_equal(ClProbX1, lProbX1)


  # Create code for a nimbleModel using dHMMo
  nc <- nimbleCode({
    x[1:5] ~ dHMMo(init[1:3], Z = Z[1:2, 1:3, 1:5],
                  T = Tt[1:3, 1:3], len = 5)

    for (i in 1:3) {
      init[i] ~ dunif(0,1)

      for (j in 1:3) {
        Tt[i,j] ~ dunif(0,1)
      }

      for (k in 1:5) {
        Z[1,i,k] ~ dunif(0,1)
        Z[2,i,k] <- 1 - Z[1,i,k]
      }
    }
  })
  # Build a nimbleModel
  m <- nimbleModel(nc, data = list(x = x1),
                   inits = list(init = init, Z = Z,
                                Tt = Tt))

  # Use the nimbleModel to calculate probabilities and compare
  m$calculate()
  MlProbX <- m$getLogProb("x")
  expect_equal(MlProbX, lProbX1)

  # Compiled model
  cm <- compileNimble(m)
  cm$calculate()
  CMlProbX <- cm$getLogProb("x")
  expect_equal(CMlProbX, lProbX1)

# Test simulation code
  set.seed(1)
  nSim <- 10
  xSim <- array(NA, dim = c(nSim, length(x1)))
  for(i in 1:nSim)
    xSim[i,] <- rHMMo(1, init, Z, Tt, len = length(x1))
  set.seed(1)
  CrHMMo <- compileNimble(rHMMo)
  CxSim <- array(NA, dim = c(nSim, length(x1)))
  for(i in 1:nSim)
    CxSim[i,] <- CrHMMo(1, init, Z, Tt, len = length(x1))
  expect_identical(xSim, CxSim)

  simNodes <- m$getDependencies(c('init', 'Z', 'Tt'), self = FALSE)
  mxSim <- array(NA, dim = c(nSim, length(x1)))
  set.seed(1)
  for(i in 1:nSim) {
    m$simulate(simNodes, includeData = TRUE)
    mxSim[i,] <- m$x
  }
  expect_identical(mxSim, xSim)

  CmxSim <- array(NA, dim = c(nSim, length(x1)))
  set.seed(1)
  for(i in 1:nSim) {
    cm$simulate(simNodes, includeData = TRUE)
    CmxSim[i,] <- cm$x
  }
  expect_identical(CmxSim, mxSim)

# # Test imputing value for all NAs
#   xNA <- rep(NA, length(x))
#   mNA <- nimbleModel(nc, data = list(x = xNA),
#                      inits = list(init = init, Z = Z,
#                                   Tt = Tt))
#   mNAConf <- configureMCMC(mNA)
#   mNAConf$addMonitors('x')
#   mNA_MCMC <- buildMCMC(mNAConf)
#   cmNA <- compileNimble(mNA, mNA_MCMC)
#   set.seed(0)
#   cmNA$mNA_MCMC$run(10)
# # Did the imputed values come back?
#   expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x[1]"])))

})

test_that("dHMM and dHMMo compatibility", {
  # len: length of the data
  len <- 4
  # Two different data samples
  x1 <- c(1, 1, 1, 2)
  x2 <- c(1, 1, 1, 1)
  # length(init) == s and sum(init) == 1; initial state probabilities
  init <- c(0.4, 0.2, 0.4)

  Z1 <- t(array(
         c(1, 0.2, 1,
           0, 0.8, 0),
         c(3, 2)))


  Z2 <- array(
         c(1, 0, 0.2, 0.8, 1, 0,
           1, 0, 0.2, 0.8, 1, 0,
           1, 0, 0.2, 0.8, 1, 0,
           1, 0, 0.2, 0.8, 1, 0),
         c(2, 3, 4))


  # Tt is time-indexed transition probabilities, s x s x t.
  Tt <- t(array(
          c(0.6, 0.3, 0.1,
            0, 0.7, 0.3,
            0, 0, 1),
          c(3,3)))

  lprob <- dHMM(x1, init, Z1, Tt, len)
  lprob_o <- dHMMo(x1, init, Z2, Tt, len)
  expect_equal(lprob, lprob_o)
})

# -----------------------------------------------------------------------------
# 3. Test that dHMM errors when input assumptions are violated
test_that("dHMM errors where expected", {

# Start with good inputs and break it one by one
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
  # T isn't square:
  badT <- t(matrix(
          c(0.6, 0.3, 0.1,
            0, 0.7, 0.3),
          ncol = length(init)))
  # Z doesn't match T:
  badZ <- t(matrix(
          c(0.6, 0.3, 0.1,
            0, 0.7, 0.3),
          ncol = length(init)))
  # Inits don't sum to 0:
  badInit <- c(0.1, 0.2, 0.2)

# Check for errors when dHMM is given bad input
  # len != length of x:
  probX <- expect_error(
              dHMM(x = x, init = init,
                   Z = Z, T = Tt,
                   len = 4, log = F))
  # T isn't square:
  probX <- expect_error(
              dHMM(x = x, init = init,
                   Z = Z, T = badT,
                   len = len, log = F))
  # Z doesn't match T:
  probX <- expect_error(
              dHMM(x = x, init = init,
                   Z = badZ, T = Tt,
                   len = len, log = F))
  # Inits don't sum to 1:
  probX <- expect_error(
              dHMM(x = x, init = badInit,
                   Z = Z, T = Tt,
                   len = len, log = F))

})

# -----------------------------------------------------------------------------
# 4. Test that dHMMo errors when input assumptions are violated

test_that("dHMMo errors where expected", {
  len <- 5
  x <- c(1, 1, 1, 2, 1)
  init <- c(0.4, 0.2, 0.4)

  Z <- array(
         c(1, 0, 0.2, 0.8, 1, 0,
           0.9, 0.1, 0.2, 0.8, 1, 0,
           1, 0, 0.2, 0.8, 1, 0,
           1, 0, 0.2, 0.8, 1, 0,
           0.95, 0.05, 0.2, 0.8, 0.5, 0.5),
         c(2, 3, 5))

  Tt <- t(matrix(
          c(0.6, 0.3, 0.1,
            0, 0.7, 0.3,
            0, 0, 1),
          ncol = length(init)))
  badT <- t(matrix(
          c(0.6, 0.3, 0.1,
            0, 0.7, 0.3),
          ncol = length(init)))
  Z_unmatched <- array(
         c(0.2, 0.8, 1, 0,
           0.2, 0.8, 1, 0,
           0.2, 0.8, 1, 0,
           0.2, 0.8, 1, 0,
           0.2, 0.8, 0.5, 0.5),
         c(2, 2, 5))
  Z_badtime <- array(
         c(1, 0, 0.2, 0.8, 1, 0,
           0.9, 0.1, 0.2, 0.8, 1, 0,
           1, 0, 0.2, 0.8, 1, 0,
           1, 0, 0.2, 0.8, 1, 0),
         c(2, 3, 4))
  badInits <- c(0.1, 0.1, 0.1)


# dHMMo tests
  # len != length of x:
  probX <- expect_error(
              dHMMo(x = x, init = init,
                   Z = Z, T = Tt,
                   len = 4, log = F))
  # T is not square
  probX <- expect_error(
              dHMMo(x = x, init = init,
                   Z = Z, T = badT,
                   len = len, log = F))
  # Z time index doesn't match len
  probX <- expect_error(
              dHMMo(x = x, init = init,
                   Z = Z_badtime, T = Tt,
                   len = len, log = F))
  # Z doesn't match T
  probX <- expect_error(
              dHMMo(x = x, init = init,
                   Z = Z_unmatched, T = Tt,
                   len = len, log = F))
  # Inits don't sum to 1
  probX <- expect_error(
              dHMMo(x = x, init = init,
                   Z = Z_unmatched, T = Tt,
                   len = len, log = F))
})



