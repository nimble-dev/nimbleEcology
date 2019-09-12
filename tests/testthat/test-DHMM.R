# Test the Dynamic Hidden Markov Model distribution nimbleFunction.

# -----------------------------------------------------------------------------
# 0. Load
# Set the context for testthat
context("Testing dDHMM-related functions.")

# -----------------------------------------------------------------------------
# 1. Test dDHMM, distribution for Dynamic Hidden Markov Model
test_that("Testing dDHMM", {
  # len: length of the data
  len <- 4
  # Two different data samples
  x1 <- c(1, 1, 1, 2)
  x2 <- c(1, 1, 1, 1)
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

  # Tt is time-indexed transition probabilities, s x s x t.
  Tt <- array(
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
    stateprob1 <- pi1 * Z[x1[i],]
    stateprob2 <- pi2 * Z[x2[i],]

    sumZ1 <- sum(stateprob1)
    sumZ2 <- sum(stateprob2)

    correctProbX1 <- correctProbX1 * sumZ1
    correctProbX2 <- correctProbX2 * sumZ2

    pi1 <- (Tt[,,i] %*% asCol(stateprob1) / sumZ1)[ ,1]
    pi2 <- (Tt[,,i] %*% asCol(stateprob2) / sumZ2)[ ,1]
  }

  # Calculate probabilities of x1 and x2 using dDHMM
  probX1 <- dDHMM(x = x1, init = init,
                  Z = Z, T = Tt,
                  len = len, log = F)
  probX2 <- dDHMM(x = x2, init = init,
                  Z = Z, T = Tt,
                  len = len, log = F)

  # Compare hand-calculated and function-calculated
  expect_equal(probX1, correctProbX1)
  expect_equal(probX2, correctProbX2)

  # Repeat for log prob
  lProbX1 <- dDHMM(x = x1, init = init,
                   Z = Z, T = Tt,
                   len = len, log = TRUE)
  lProbX2 <- dDHMM(x = x2, init = init,
                   Z = Z, T = Tt,
                   len = len, log = TRUE)

  expect_equal(lProbX1, log(correctProbX1))
  expect_equal(lProbX2, log(correctProbX2))

  # Repeat for the compiled function
  CdDHMM <- compileNimble(dDHMM)
  CprobX1 <- CdDHMM(x = x1, init = init,
                    Z = Z, T = Tt,
                    len = len, log = F)
  expect_equal(CprobX1, probX1)

  ClProbX1 <- CdDHMM(x = x1, init = init,
                     Z = Z, T = Tt,
                     len = len, log = TRUE)
  expect_equal(ClProbX1, lProbX1)

  # Check if the random generator rHMM works
  # These values aren't necessarily true, but this will break if the function
  # is changed and needs to be re-tested
  set.seed(1111)
  oneSim <- rDHMM(1, init, Z, Tt, len = len)
  expect_equal(oneSim, c(2, 2, 1, 1))

    # Create code for a nimbleModel using the distribution
  nc <- nimbleCode({
    x[1:4] ~ dDHMM(init[1:3], Z = Z[1:2,1:3],
                   T = Tt[1:3, 1:3, 1:4], len = 4)

    for (i in 1:3) {
      init[i] ~ dunif(0,1)

      for (j in 1:3) {
        for (t in 1:4) {
          Tt[i,j,t] ~ dunif(0,1)
        }
      }

      Z[1,i] ~ dunif(0,1)
      Z[2,i] <- 1 - Z[1,i]
    }
  })

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

  # Check simulate
  set.seed(2468)
  cm$simulate('x')
  expect_equal(cm$x, x1)

})


# -----------------------------------------------------------------------------
# 2. Test dDHMMo, distribution for Dynamic Hidden Markov Model with time-
#    indexed observation probabilities
test_that("Testing dDHMMo", {
  # len: length of the data
  len <- 4
  # Two different data samples
  x1 <- c(1, 1, 1, 2)
  x2 <- c(1, 1, 1, 1)
  # length(init) == s and sum(init) == 1; initial state probabilities
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
           0.95, 0.05, 0.2, 0.8, 0.5, 0.5),
         c(2, 3, 4))


  # Tt is time-indexed transition probabilities, s x s x t.
  Tt <- array(
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
    stateprob1 <- pi1 * Z[x1[i],,i]
    stateprob2 <- pi2 * Z[x2[i],,i]

    sumZ1 <- sum(stateprob1)
    sumZ2 <- sum(stateprob2)

    correctProbX1 <- correctProbX1 * sumZ1
    correctProbX2 <- correctProbX2 * sumZ2

    pi1 <- (Tt[,,i] %*% asCol(stateprob1) / sumZ1)[ ,1]
    pi2 <- (Tt[,,i] %*% asCol(stateprob2) / sumZ2)[ ,1]
  }

  # Calculate probabilities of x1 and x2 using dDHMM
  probX1 <- dDHMMo(x = x1, init = init,
                   Z = Z, T = Tt,
                   len = len, log = F)
  probX2 <- dDHMMo(x = x2, init = init,
                   Z = Z, T = Tt,
                   len = len, log = F)

  # Compare hand-calculated and function-calculated
  expect_equal(probX1, correctProbX1)
  expect_equal(probX2, correctProbX2)

  # Repeat for log prob
  lProbX1 <- dDHMMo(x = x1, init = init,
                    Z = Z, T = Tt,
                    len = len, log = TRUE)
  lProbX2 <- dDHMMo(x = x2, init = init,
                    Z = Z, T = Tt,
                    len = len, log = TRUE)

  expect_equal(lProbX1, log(correctProbX1))
  expect_equal(lProbX2, log(correctProbX2))

  # Repeat for the compiled function
  CdDHMMo <- compileNimble(dDHMMo)
  CprobX1 <- CdDHMMo(x = x1, init = init,
                     Z = Z, T = Tt,
                     len = len, log = F)
  expect_equal(CprobX1, probX1)

  ClProbX1 <- CdDHMMo(x = x1, init = init,
                      Z = Z, T = Tt,
                      len = len, log = TRUE)
  expect_equal(ClProbX1, lProbX1)

  # Check if the random generator rHMM works
  # These values aren't necessarily true, but this will break if the function
  # is changed and needs to be re-tested
  set.seed(1111)
  oneSim <- rDHMMo(1, init, Z, Tt, len = len)
  expect_equal(oneSim, c(2, 2, 1, 1))

  # Create code for a nimbleModel using the distribution
  nc <- nimbleCode({
    x[1:4] ~ dDHMMo(init[1:3], Z = Z[1:2,1:3,1:4],
                    T = Tt[1:3, 1:3, 1:4], len = 4)

    for (i in 1:3) {
      init[i] ~ dunif(0,1)

      for (t in 1:4) {
        Z[1,i,t] ~ dunif(0,1)
        Z[2,i,t] <- 1 - Z[1,i,t]
        for (j in 1:3) {
          Tt[i,j,t] ~ dunif(0,1)
        }
      }
    }
  })

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

  # Check simulate
  set.seed(2468)
  cm$simulate('x')
  expect_equal(cm$x, x1)

})

# TODO:
# Check that time-indexed with all the same probs for each time matches
# non-time-indexed

test_that("dDHMM and dDHMMo compatibility", {
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
  Tt <- array(
          c(0.6, 0, 0, 0.3, 0.7, 0.25, 0.1, 0.3, 0.75,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0, 0.2, 0.7, 0, 0.2, 0.3, 1,
            0.6, 0, 0.2, 0.3, 0.7, 0, 0.1, 0.3, 0.8
            ),
          c(3,3,4))

  lprob <- dDHMM(x1, init, Z1, Tt, len)
  lprob_o <- dDHMMo(x1, init, Z2, Tt, len)
  expect_equal(lprob, lprob_o)
})
