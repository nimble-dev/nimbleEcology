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
Z <- t(matrix(
       c(1, 0.2, 1,
         0, 0.8, 0),
       nrow = length(init))) # I'm doing the transpose to specify by rows, but remember that nxxx is backwards

# Tt is transition probabilities, s x s.
Tt <- t(matrix(
        c(0.6, 0.3, 0.1,
          0, 0.7, 0.3,
          0, 0, 1),
        ncol = length(init)))


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


})



