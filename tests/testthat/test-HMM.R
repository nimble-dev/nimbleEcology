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

  # probObs is observation probabilities, dim o x s where o is possible observation
  # states (domain of x) and s is possible true states. probObs says, for each true
  # system state (row), what is the corresponding probability of observing each
  # response (col)
  probObs <- t(array(
         c(1, 0,
           0, 1,
           0.8, 0.2),
         c(2, 3)))

  # probTrans is transition probabilities, s x s.
  probTrans <- t(array(
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
    stateprob1 <- pi1 * probObs[, x1[i]]
    stateprob2 <- pi2 * probObs[, x2[i]]

    sumZ1 <- sum(stateprob1)
    sumZ2 <- sum(stateprob2)

    correctProbX1 <- correctProbX1 * sumZ1
    correctProbX2 <- correctProbX2 * sumZ2

    pi1 <- (asRow(stateprob1) %*% probTrans[, ]/sumZ1)[1, ]
    pi2 <- (asRow(stateprob2) %*% probTrans[, ]/sumZ2)[1, ]
  }

  # Calculate probabilities of x1 and x2 using dHMM
  probX1 <- dHMM(x = x1, init = init,
                 probObs = probObs, probTrans = probTrans,
                 len = len, log = F)
  probX2 <- dHMM(x = x2, init = init,
                 probObs = probObs, probTrans = probTrans,
                 len = len, log = F)

  # Compare hand-calculated and function-calculated
  expect_equal(probX1, correctProbX1)
  expect_equal(probX2, correctProbX2)

  # Repeat for log prob
  lProbX1 <- dHMM(x = x1, init = init,
                  probObs = probObs, probTrans = probTrans,
                  len = len, log = TRUE)
  lProbX2 <- dHMM(x = x2, init = init,
                  probObs = probObs, probTrans = probTrans,
                  len = len, log = TRUE)

  expect_equal(lProbX1, log(correctProbX1))
  expect_equal(lProbX2, log(correctProbX2))

  # Repeat for the compiled function
  CdHMM <- compileNimble(dHMM)
  CprobX1 <- CdHMM(x = x1, init = init,
                   probObs = probObs, probTrans = probTrans,
                   len = len, log = FALSE)
  expect_equal(CprobX1, probX1)

  ClProbX1 <- CdHMM(x = x1, init = init,
                    probObs = probObs, probTrans = probTrans,
                    len = len, log = TRUE)
  expect_equal(ClProbX1, lProbX1)

  # Create code for a nimbleModel using the distribution
  nc <- nimbleCode({
    x[1:5] ~ dHMM(init[1:3], probObs = probObs[1:3,1:2],
                  probTrans = probTrans[1:3, 1:3], len = 5)
  })

  # Create a nimbleModel using the distribution
  m <- nimbleModel(nc, data = list(x = x1),
                   inits = list(init = init, probObs = probObs,
                                probTrans = probTrans))
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
    xSim[i,] <- rHMM(1, init, probObs, probTrans, len = length(x1))
  set.seed(1)
  CrHMM <- compileNimble(rHMM)
  CxSim <- array(NA, dim = c(nSim, length(x1)))
  for(i in 1:nSim)
    CxSim[i,] <- CrHMM(1, init, probObs, probTrans, len = length(x1))
  expect_identical(xSim, CxSim)

  simNodes <- m$getDependencies(c('init', 'probObs', 'probTrans'), self = FALSE)
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

#   # Test imputing value for all NAs
#   xNA <- rep(NA, length(x))
#   mNA <- nimbleModel(nc, data = list(x = xNA),
#                      inits = list(init = init, probObs = probObs,
#                                   probTrans = probTrans))
#   mNAConf <- configureMCMC(mNA)
#   mNAConf$addMonitors('x')
#   mNA_MCMC <- buildMCMC(mNAConf)
#   cmNA <- compileNimble(mNA, mNA_MCMC)
#   set.seed(0)
#   cmNA$mNA_MCMC$run(10)
# # Did the imputed values come back?
#   expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x[1]"])))
})


# -----------------------------------------------------------------------------
# 2. Test dHMMo, distribution for Hidden Markov Model
test_that("dHMMo works", {
  len <- 5
  x1 <- c(1, 1, 1, 2, 1)
  x2 <- c(1, 1, 1, 1, 1)
  # length(init) == s
  init <- c(0.4, 0.2, 0.4)

  # probObs is observation probabilities, dim o x s x t where o is possible
  # observation states (domain of x), s is possible true system states, and t is
  # time index. probObs says, for each true system state (row), what is the
  # corresponding probability of observing each response (col) at time t (3rd
  # dim)
  probObs <- array(
         c(1, 0, 0.8, 0, 1, 0.2,
           1, 0, 0.8, 0, 1, 0.2,
           0.9, 0, 0.8, 0.1, 1, 0.2,
           1, 0.1, 0.8, 0, 0.9, 0.2,
           1, 0, 0.7, 0, 1, 0.3),
         c(3, 2, 5))

  # probTrans is transition probabilities, s x s.
  probTrans <- t(array(
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
    stateprob1 <- pi1 * probObs[,x1[i],i]
    stateprob2 <- pi2 * probObs[,x2[i],i]

    sumZ1 <- sum(stateprob1)
    sumZ2 <- sum(stateprob2)

    correctProbX1 <- correctProbX1 * sumZ1
    correctProbX2 <- correctProbX2 * sumZ2

    pi1 <- (asRow(stateprob1) %*% probTrans[, ]/sumZ1)[1, ]
    pi2 <- (asRow(stateprob2) %*% probTrans[, ]/sumZ2)[1, ]
  }

  # Calculate probabilities of xs using the function dHMMo
  probX1 <- dHMMo(x = x1, init = init,
                  probObs = probObs, probTrans = probTrans,
                  len = len, log = F)
  probX2 <- dHMMo(x = x2, init = init,
                  probObs = probObs, probTrans = probTrans,
                  len = len, log = F)

  # Compare function and manual
  expect_equal(probX1, correctProbX1)
  expect_equal(probX2, correctProbX2)

  # Repeat for log prob
  lProbX1 <- dHMMo(x = x1, init = init,
                   probObs = probObs, probTrans = probTrans,
                   len = len, log = TRUE)
  lProbX2 <- dHMMo(x = x2, init = init,
                   probObs = probObs, probTrans = probTrans,
                   len = len, log = TRUE)

  expect_equal(lProbX1, log(correctProbX1))
  expect_equal(lProbX2, log(correctProbX2))

  # Repeat for compiled nimbleFunction
  CdHMMo <- compileNimble(dHMMo)
  CprobX1 <- CdHMMo(x = x1, init = init,
                    probObs = probObs, probTrans = probTrans,
                    len = len, log = FALSE)
  expect_equal(CprobX1, probX1)

  ClProbX1 <- CdHMMo(x = x1, init = init,
                     probObs = probObs, probTrans = probTrans,
                     len = len, log = TRUE)
  expect_equal(ClProbX1, lProbX1)


  # Create code for a nimbleModel using dHMMo
  nc <- nimbleCode({
    x[1:5] ~ dHMMo(init[1:3], probObs = probObs[1:3, 1:2, 1:5],
                  probTrans = probTrans[1:3, 1:3], len = 5)
  })
  # Build a nimbleModel
  m <- nimbleModel(nc, data = list(x = x1),
                   inits = list(init = init, probObs = probObs,
                                probTrans = probTrans))

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
    xSim[i,] <- rHMMo(1, init, probObs, probTrans, len = length(x1))
  set.seed(1)
  CrHMMo <- compileNimble(rHMMo)
  CxSim <- array(NA, dim = c(nSim, length(x1)))
  for(i in 1:nSim)
    CxSim[i,] <- CrHMMo(1, init, probObs, probTrans, len = length(x1))
  expect_identical(xSim, CxSim)

  simNodes <- m$getDependencies(c('init', 'probObs', 'probTrans'), self = FALSE)
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
#                      inits = list(init = init, probObs = probObs,
#                                   probTrans = probTrans))
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
  len <- 5
  # Two different data samples
  x1 <- c(1, 1, 1, 2, 1)
  x2 <- c(1, 1, 1, 1, 1)
  # length(init) == s and sum(init) == 1; initial state probabilities
  init <- c(0.4, 0.2, 0.4)

  probObs1 <- t(array(
         c(1, 0,
           0, 1,
           0.8, 0.2),
         c(2, 3)))


  probObs2 <- array(
         c(1, 0, 0.8, 0, 1, 0.2,
           1, 0, 0.8, 0, 1, 0.2,
           1, 0, 0.8, 0, 1, 0.2,
           1, 0, 0.8, 0, 1, 0.2,
           1, 0, 0.8, 0, 1, 0.2),
         c(3, 2, 5))

  # probTrans is time-indexed transition probabilities, s x s x t.
  probTrans <- t(array(
          c(0.6, 0.3, 0.1,
            0, 0.7, 0.3,
            0, 0, 1),
          c(3,3)))

  lprob <- dHMM(x1, init, probObs1, probTrans, len)
  lprob_o <- dHMMo(x1, init, probObs2, probTrans, len)
  expect_equal(lprob, lprob_o)
})

# -----------------------------------------------------------------------------
# 3. Test that dHMM errors when input assumptions are violated
test_that("dHMM errors where expected", {

# Start with good inputs and break it one by one
  len <- 5
  x <- c(1, 1, 1, 2, 1)
  init <- c(0.4, 0.2, 0.4)

  badprobObs <- t(matrix(
         c(1, 0.2, 1,
           0, 0.8, 0),
         nrow = length(init)))

  probTrans <- t(matrix(
          c(0.6, 0.3, 0.1,
            0, 0.7, 0.3,
            0, 0, 1),
          ncol = length(init)))
  # T isn't square:
  badT <- t(matrix(
          c(0.6, 0.3, 0.1,
            0, 0.7, 0.3),
          ncol = length(init)))
  # probObs doesn't match T:
  probObs <- t(matrix(
          c(0.6, 0.3, 0.1,
            0, 0.7, 0.3),
          ncol = length(init)))
  # Inits don't sum to 0:
  badInit <- c(0.1, 0.2, 0.2)

# Check for errors when dHMM is given bad input
  # len != length of x:
  probX <- expect_error(
              dHMM(x = x, init = init,
                   probObs = probObs, probTrans = probTrans,
                   len = 4, log = F))
  # T isn't square:
  probX <- expect_error(
              dHMM(x = x, init = init,
                   probObs = probObs, probTrans =badT,
                   len = len, log = F))
  # probObs doesn't match T:
  probX <- expect_error(
              dHMM(x = x, init = init,
                   probObs = badprobObs, probTrans = probTrans,
                   len = len, log = F))
  # Inits don't sum to 1:
  probX <- expect_error(
              dHMM(x = x, init = badInit,
                   probObs = probObs, probTrans = probTrans,
                   len = len, log = F))

})

# -----------------------------------------------------------------------------
# 4. Test that dHMMo errors when input assumptions are violated

test_that("dHMMo errors where expected", {
  len <- 5
  x <- c(1, 1, 1, 2, 1)
  init <- c(0.4, 0.2, 0.4)

  probObs <- array(
         c(1, 0, 0.8, 0, 1, 0.2,
           1, 0, 0.8, 0, 1, 0.2,
           1, 0, 0.8, 0, 1, 0.2,
           1, 0, 0.8, 0, 1, 0.2,
           1, 0, 0.8, 0, 1, 0.2),
         c(3, 2, 5))

  probTrans <- t(matrix(
          c(0.6, 0.3, 0.1,
            0, 0.7, 0.3,
            0, 0, 1),
          ncol = length(init)))
  badT <- t(matrix(
          c(0.6, 0.3, 0.1,
            0, 0.7, 0.3),
          ncol = length(init)))
  probObs_unmatched <- array(
         c(0.2, 0.8, 1, 0,
           0.2, 0.8, 1, 0,
           0.2, 0.8, 1, 0,
           0.2, 0.8, 1, 0,
           0.2, 0.8, 0.5, 0.5),
         c(2, 2, 5))
  probObs_badtime <- array(
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
                   probObs = probObs, probTrans = probTrans,
                   len = 4, log = F))
  # T is not square
  probX <- expect_error(
              dHMMo(x = x, init = init,
                   probObs = probObs, probTrans =badT,
                   len = len, log = F))
  # probObs time index doesn't match len
  probX <- expect_error(
              dHMMo(x = x, init = init,
                   probObs = probObs_badtime, probTrans = probTrans,
                   len = len, log = F))
  # probObs doesn't match T
  probX <- expect_error(
              dHMMo(x = x, init = init,
                   probObs = probObs_unmatched, probTrans = probTrans,
                   len = len, log = F))
  # Inits don't sum to 1
  probX <- expect_error(
              dHMMo(x = x, init = init,
                   probObs = probObs_unmatched, probTrans = probTrans,
                   len = len, log = F))
})



