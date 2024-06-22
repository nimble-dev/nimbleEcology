# Test the Dynamic Hidden Markov Model distribution nimbleFunction.

# -----------------------------------------------------------------------------
# 0. Load
# -----------------------------------------------------------------------------
# 1. Test dDHMM, distribution for Dynamic Hidden Markov Model
test_that("Testing dDHMM", {
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

  # probTrans is time-indexed transition probabilities, s x s x t.
  probTrans <- array(
          c(0.6, 0, 0, 0.3, 0.7, 0.25, 0.1, 0.3, 0.75,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0.2, 0.3, 0.7, 0, 0.1, 0.3, 0.8
            ),
          c(3,3,4))

  correctProbX1 <- 1
  correctProbX2 <- 1
  pi1 <- init
  pi2 <- init
  for (i in 1:len) {
    stateprob1 <- pi1 * probObs[,x1[i]]
    stateprob2 <- pi2 * probObs[,x2[i]]

    sumZ1 <- sum(stateprob1)
    sumZ2 <- sum(stateprob2)

    correctProbX1 <- correctProbX1 * sumZ1
    correctProbX2 <- correctProbX2 * sumZ2

    if (i != len) pi1 <- (asRow(stateprob1) %*% probTrans[,,i]/sumZ1)[1, ]
    if (i != len) pi2 <- (asRow(stateprob2) %*% probTrans[,,i]/sumZ2)[1, ]
  }

  # Calculate probabilities of x1 and x2 using dDHMM
  probX1 <- dDHMM(x = x1, init = init,
                  probObs = probObs, probTrans = probTrans,
                  len = len, log = F)
  probX2 <- dDHMM(x = x2, init = init,
                  probObs = probObs, probTrans = probTrans,
                  len = len, log = F)

  # Compare hand-calculated and function-calculated
  expect_equal(probX1, correctProbX1)
  expect_equal(probX2, correctProbX2)

  # Repeat for log prob
  lProbX1 <- dDHMM(x = x1, init = init,
                   probObs = probObs, probTrans = probTrans,
                   len = len, log = TRUE)
  lProbX2 <- dDHMM(x = x2, init = init,
                   probObs = probObs, probTrans = probTrans,
                   len = len, log = TRUE)

  expect_equal(lProbX1, log(correctProbX1))
  expect_equal(lProbX2, log(correctProbX2))

  call_dDHMM <- nimbleFunction(
    run = function(x=double(1), init=double(1), probObs=double(2),
                   probTrans=double(3), len=integer(0,default=0),
                   checkRowSums = integer(0,default=1),
                   log=integer(0, default=0)) {
      return(dDHMM(x,init,probObs,probTrans,len,checkRowSums,log))
      returnType(double())
    })

  # Repeat for the compiled function
  CdDHMM <- compileNimble(call_dDHMM)
  CprobX1 <- CdDHMM(x = x1, init = init,
                    probObs = probObs, probTrans = probTrans,
                    len = len, log = F)
  expect_equal(CprobX1, probX1)

  ClProbX1 <- CdDHMM(x = x1, init = init,
                     probObs = probObs, probTrans = probTrans,
                     len = len, log = TRUE)
  expect_equal(ClProbX1, lProbX1)

    # Create code for a nimbleModel using the distribution
  nc <- nimbleCode({
    x[1:5] ~ dDHMM(init[1:3], probObs = probObs[1:3,1:2],
                   probTrans = probTrans[1:3, 1:3, 1:4], len = 5, checkRowSums = 1)
  })

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
    xSim[i,] <- rDHMM(1, init, probObs, probTrans, len = length(x1))
  set.seed(1)
  CrDHMM <- compileNimble(rDHMM)
  CxSim <- array(NA, dim = c(nSim, length(x1)))
  for(i in 1:nSim)
    CxSim[i,] <- CrDHMM(1, init, probObs, probTrans, len = length(x1))
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


# -----------------------------------------------------------------------------
# 2. Test dDHMMo, distribution for Dynamic Hidden Markov Model with time-
#    indexed observation probabilities
test_that("Testing dDHMMo", {
  # len: length of the data
  len <- 5
  # Two different data samples
  x1 <- c(1, 1, 1, 2, 1)
  x2 <- c(1, 1, 1, 1, 1)
  # length(init) == s and sum(init) == 1; initial state probabilities
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


  # probTrans is time-indexed transition probabilities, s x s x t.
  probTrans <- array(
          c(0.6, 0, 0, 0.3, 0.7, 0.25, 0.1, 0.3, 0.75,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0.2, 0.3, 0.7, 0, 0.1, 0.3, 0.8
            ),
          c(3,3,4))

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

    if (i != len) pi1 <- (asRow(stateprob1) %*% probTrans[,,i]/sumZ1)[1, ]
    if (i != len) pi2 <- (asRow(stateprob2) %*% probTrans[,,i]/sumZ2)[1, ]
  }

  # Calculate probabilities of x1 and x2 using dDHMM
  probX1 <- dDHMMo(x = x1, init = init,
                   probObs = probObs, probTrans = probTrans,
                   len = len, log = F)
  probX2 <- dDHMMo(x = x2, init = init,
                   probObs = probObs, probTrans = probTrans,
                   len = len, log = F)

  # Compare hand-calculated and function-calculated
  expect_equal(probX1, correctProbX1)
  expect_equal(probX2, correctProbX2)

  # Repeat for log prob
  lProbX1 <- dDHMMo(x = x1, init = init,
                    probObs = probObs, probTrans = probTrans,
                    len = len, log = TRUE)
  lProbX2 <- dDHMMo(x = x2, init = init,
                    probObs = probObs, probTrans = probTrans,
                    len = len, log = TRUE)

  expect_equal(lProbX1, log(correctProbX1))
  expect_equal(lProbX2, log(correctProbX2))

  # Repeat for the compiled function
  call_dDHMMo <- nimbleFunction(
    run = function(x=double(1), init=double(1), probObs=double(3),
                   probTrans=double(3), len=integer(0,default=0),
                   checkRowSums = integer(0,default=1),
                   log=integer(0, default=0)) {
      return(dDHMMo(x,init,probObs,probTrans,len,checkRowSums,log))
      returnType(double())
    })

  CdDHMMo <- compileNimble(call_dDHMMo)
  CprobX1 <- CdDHMMo(x = x1, init = init,
                     probObs = probObs, probTrans = probTrans,
                     len = len, log = F)
  expect_equal(CprobX1, probX1)

  ClProbX1 <- CdDHMMo(x = x1, init = init,
                      probObs = probObs, probTrans = probTrans,
                      len = len, log = TRUE)
  expect_equal(ClProbX1, lProbX1)

  # Create code for a nimbleModel using the distribution
  nc <- nimbleCode({
    x[1:5] ~ dDHMMo(init[1:3], probObs = probObs[1:3,1:2,1:5],
                    probTrans = probTrans[1:3, 1:3, 1:4], len = 5, checkRowSums = 1)
  })

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
    xSim[i,] <- rDHMMo(1, init, probObs, probTrans, len = length(x1))
  set.seed(1)
  CrDHMMo <- compileNimble(rDHMMo)
  CxSim <- array(NA, dim = c(nSim, length(x1)))
  for(i in 1:nSim)
    CxSim[i,] <- CrDHMMo(1, init, probObs, probTrans, len = length(x1))
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

test_that("dDHMM and dDHMMo compatibility", {
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
  probTrans <- array(
          c(0.6, 0, 0, 0.3, 0.7, 0.25, 0.1, 0.3, 0.75,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0.2, 0.3, 0.7, 0, 0.1, 0.3, 0.8
            ),
          c(3,3,4))

  lprob <- dDHMM(x1, init, probObs1, probTrans, len, log = TRUE)
  lprob_o <- dDHMMo(x1, init, probObs2, probTrans, len, log = TRUE)
  expect_equal(lprob, lprob_o)
})


test_that("dDHMM errors where expected", {
  len <- 5
  x <- c(1, 1, 1, 2, 1)
  init <- c(0.4, 0.2, 0.4)

  probObs <- t(array(
         c(1, 0,
           0, 1,
           0.8, 0.2),
         c(2, 3)))


  probTrans <- array(
          c(0.6, 0, 0, 0.3, 0.7, 0.25, 0.1, 0.3, 0.75,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0.2, 0.3, 0.7, 0, 0.1, 0.3, 0.8
            ),
          c(3,3,4))
  badT <- array(
          c(0.6, 0, 0, 0.3, 0.7, 0.25, 0.1, 0.3, 0.75,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0.2, 0.3, 0.7, 0, 0.1, 0.3, 0.8
            ),
          c(10,3,4))
  probObs_unmatched <- t(array(
         c(1, 0,
           0, 1,
           0.8, 0.2),
         c(1, 6)))

  badInits <- c(0.1, 0.1, 0.1)


# dHMMo tests
  # len != length of x:
  probX <- expect_error(
              dDHMM(x = x, init = init,
                   probObs = probObs, probTrans = probTrans,
                   len = 4, log = F))
  # T is not square
  probX <- expect_error(
              dDHMM(x = x, init = init,
                   probObs = probObs, probTrans =badT,
                   len = len, log = F))
  # probObs doesn't match T
  probX <- expect_error(
              dDHMM(x = x, init = init,
                   probObs = probObs_unmatched, probTrans = probTrans,
                   len = len, log = F))
  # Inits don't sum to 1
  probX <- expect_error(
              dDHMM(x = x, init = badInits,
                   probObs = probObs, probTrans = probTrans,
                   len = len, log = F))

  # Bad sums for probObs:
  bpo2 <- probObs
  bpo2[1,] <- 0
  probX <- expect_error(
            dDHMM(x = x, init = init,
                 probObs = bpo2, probTrans = probTrans,
                 len = len, log = F))

  # Bad sums for probTrans:
  bpt2 <- probTrans
  bpt2[1,,1] <- 0
  probX <- expect_error(
            dDHMM(x = x, init = init,
                 probObs = probObs, probTrans = bpt2,
                 len = len, log = F))
})


test_that("dDHMMo errors where expected", {
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

  probTrans <- array(
          c(0.6, 0, 0, 0.3, 0.7, 0.25, 0.1, 0.3, 0.75,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0.2, 0.3, 0.7, 0, 0.1, 0.3, 0.8
            ),
          c(3,3,4))
  badT <- array(
          c(0.6, 0, 0, 0.3, 0.7, 0.25, 0.1, 0.3, 0.75,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0.2, 0.3, 0.7, 0, 0.1, 0.3, 0.8
            ),
          c(1,9,4))
  badT2 <- array(
          c(0.8, 0, 0, 0.3, 0.7, 0.25, 0.1, 0.3, 0.75,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0.2, 0.3, 0.7, 0, 0.1, 0.3, 0.8
            ),
          c(3,3,4))

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
              dDHMMo(x = x, init = init,
                   probObs = probObs, probTrans = probTrans,
                   len = 4, log = F))
  # T is not square
  probX <- expect_error(
              dDHMMo(x = x, init = init,
                   probObs = probObs, probTrans =badT,
                   len = len, log = F))
  # probObs time index doesn't match len
  probX <- expect_error(
              dDHMMo(x = x, init = init,
                   probObs = probObs_badtime, probTrans = probTrans,
                   len = len, log = F))
  # probObs doesn't match T
  probX <- expect_error(
              dDHMMo(x = x, init = init,
                   probObs = probObs_unmatched, probTrans = probTrans,
                   len = len, log = F))
  # Inits don't sum to 1
  probX <- expect_error(
              dDHMMo(x = x, init = badInit,
                   probObs = probObs, probTrans = probTrans,
                   len = len, log = F))

  # Bad sums for probObs:
  bpo2 <- probObs
  bpo2[1,,] <- 0
  probX <- expect_error(
            dDHMMo(x = x, init = init,
                 probObs = bpo2, probTrans = probTrans,
                 len = len, log = F))

  # Bad sums for probTrans:
  probX <- expect_error(
            dDHMMo(x = x, init = init,
                 probObs = probObs, probTrans = badT2,
                 len = len, log = F))
})





test_that("rDHMM errors where expected", {
  len <- 5
  x <- c(1, 1, 1, 2, 1)
  init <- c(0.4, 0.2, 0.4)

  probObs <- t(array(
         c(1, 0,
           0, 1,
           0.8, 0.2),
         c(2, 3)))


  probTrans <- array(
          c(0.6, 0, 0, 0.3, 0.7, 0.25, 0.1, 0.3, 0.75,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0.2, 0.3, 0.7, 0, 0.1, 0.3, 0.8
            ),
          c(3,3,4))
  badT <- array(
          c(0.6, 0, 0, 0.3, 0.7, 0.25, 0.1, 0.3, 0.75,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0.2, 0.3, 0.7, 0, 0.1, 0.3, 0.8
            ),
          c(10,3,4))
  probObs_unmatched <- t(array(
         c(1, 0,
           0, 1,
           0.8, 0.2),
         c(1, 6)))

  badInits <- c(0.1, 0.1, 0.1)


# dHMMo tests
  # T is not square
  probX <- expect_error(
              rDHMM(init = init,
                   probObs = probObs, probTrans = badT,
                   len = len))
  # probObs doesn't match T
  probX <- expect_error(
              rDHMM(init = init,
                   probObs = probObs_unmatched, probTrans = probTrans,
                   len = len))
  # Inits don't sum to 1
  probX <- expect_error(
              rDHMM(init = badInits,
                   probObs = probObs, probTrans = probTrans,
                   len = len))

  # Bad sums for probObs:
  bpo2 <- probObs
  bpo2[1,] <- 0
  probX <- expect_error(
            rDHMM(init = init,
                 probObs = bpo2, probTrans = probTrans,
                 len = len))

  # Bad sums for probTrans:
  bpt2 <- probTrans
  bpt2[1,,1] <- 0
  probX <- expect_error(
            rDHMM(init = init,
                 probObs = probObs, probTrans = bpt2,
                 len = len))
})


test_that("rDHMMo errors where expected", {
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

  probTrans <- array(
          c(0.6, 0, 0, 0.3, 0.7, 0.25, 0.1, 0.3, 0.75,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0.2, 0.3, 0.7, 0, 0.1, 0.3, 0.8
            ),
          c(3,3,4))
  badT <- array(
          c(0.6, 0, 0, 0.3, 0.7, 0.25, 0.1, 0.3, 0.75,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0.2, 0.3, 0.7, 0, 0.1, 0.3, 0.8
            ),
          c(1,9,4))
  badT2 <- array(
          c(0.8, 0, 0, 0.3, 0.7, 0.25, 0.1, 0.3, 0.75,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0.2, 0.3, 0.7, 0, 0.1, 0.3, 0.8
            ),
          c(3,3,4))

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


# rDHMMo tests
  # T is not square
  probX <- expect_error(
              rDHMMo(init = init,
                   probObs = probObs, probTrans =badT,
                   len = len))
  # probObs time index doesn't match len
  probX <- expect_error(
              rDHMMo(init = init,
                   probObs = probObs_badtime, probTrans = probTrans,
                   len = len))
  # probObs doesn't match T
  probX <- expect_error(
              rDHMMo(init = init,
                   probObs = probObs_unmatched, probTrans = probTrans,
                   len = len))
  # Inits don't sum to 1
  probX <- expect_error(
              rDHMMo(init = badInit,
                   probObs = probObs, probTrans = probTrans,
                   len = len))

  # Bad sums for probObs:
  bpo2 <- probObs
  bpo2[1,,] <- 0
  probX <- expect_error(
            rDHMMo(init = init,
                 probObs = bpo2, probTrans = probTrans,
                 len = len))

  # Bad sums for probTrans:
  probX <- expect_error(
            rDHMMo(init = init,
                 probObs = probObs, probTrans = badT2,
                 len = len))
})
