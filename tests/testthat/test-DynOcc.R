
library(testthat)
library(nimble)
source("R/dDynOcc.R")
context("Testing dDynOcc-related functions.")


test_that("dDynOcc_vv works", {
  x <- matrix(c(0,0,1,0,1,
                1,1,1,1,0), nrow = 2)
  nrep <- c(5, 5)
  psi1 <- 0.7
  phi <- c(0.4, 0.8)
  gamma <- c(0.2, 0.1)
  p <- matrix(rep(0.8, 10), nrow = 2)

  probX <- dDynOcc_vv(x, nrep, psi1, phi, gamma, p, log = FALSE)
  lProbX <- dDynOcc_vv(x, nrep, psi1, phi, gamma, p, log = TRUE)

  ProbOccNextTime <- psi1
  ll <- 0
  nyears <- dim(x)[1]
  if (nyears >= 1) {
    for (t in 1:nyears) {
      if (nrep[t] > 0) {
        numObs <- sum(x[t,1:nrep[t]])
        if (numObs < 0) {
          print("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
          stop("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
        }
        ProbOccAndCount <- ProbOccNextTime *
            exp(sum(dbinom(x[t,1:nrep[t]],
                           size = 1, p = p[t,1:nrep[t]], log = 1)))
        ProbUnoccAndCount <- (1-ProbOccNextTime) * (numObs == 0)
        ProbCount <- ProbOccAndCount + ProbUnoccAndCount
        ProbOccGivenCount <- ProbOccAndCount / ProbCount
        ll <- ll + log(ProbCount)
        if(t < nyears)
            ProbOccNextTime <- ProbOccGivenCount * phi[t] +
                (1-ProbOccGivenCount) * gamma[t]
      } else {
        ## If there were no observations in year t,
        ## simply propagate probability of occupancy forward
        if (t < nyears)
            ProbOccNextTime <- ProbOccNextTime * phi[t] +
                (1-ProbOccNextTime) * gamma[t]
      }
    }
  }

  lCorrectProbX <- ll

  expect_equal(exp(lCorrectProbX), probX)
  expect_equal(lCorrectProbX, lProbX)

  CdDynOcc <- compileNimble(dDynOcc_vv)
  CprobX <- CdDynOcc(x, nrep, psi1, phi, gamma, p, log = FALSE)
  expect_equal(CprobX, probX)

  ClProbX <- CdDynOcc(x, nrep, psi1, phi, gamma, p, log = TRUE)
  expect_equal(ClProbX, lProbX)



  nc <- nimbleCode({

    x[1:2, 1:5] ~ dDynOcc_vv(nrep[1:2], psi1, phi[1:2], gamma[1:2], p[1:2,1:5])

    psi1 ~ dunif(0,1)

    for (i in 1:2) {
      phi[i] ~ dunif(0,1)
      gamma[i] ~ dunif(0,1)
    }

    for (i in 1:2) {
      for (j in 1:5) {
        p[i,j] ~ dunif(0,1)
      }
    }
  })

  m <- nimbleModel(nc, data = list(x = x, nrep = nrep),
                   # constants = list(nrep = nrep),
                   inits = list(p = p, phi = phi,
                                psi1 = psi1, gamma = gamma))

  m$calculate()
  MlProbX <- m$getLogProb("x")
  expect_equal(MlProbX, lProbX)

  cm <- compileNimble(m, showCompilerOutput = TRUE)
  cm$calculate()
  CMlProbX <- cm$getLogProb("x")
  expect_equal(CMlProbX, lProbX)

  set.seed(2468)
  cm$simulate('x')
  expect_equal(cm$x, matrix(c(0, 0, 1, 0, 1, 1, 1, 1, 1, 0), nrow = 2))
})


test_that("dDynOcc_vs works", {
  x <- matrix(c(0,0,1,0,1,
                1,1,1,1,0), nrow = 2)
  nrep <- c(5, 5)
  psi1 <- 0.7
  phi <- c(0.4, 0.8)
  gamma <- 0.5
  p <- matrix(rep(0.8, 10), nrow = 2)

  probX <- dDynOcc_vs(x, nrep, psi1, phi, gamma, p, log = FALSE)
  lProbX <- dDynOcc_vs(x, nrep, psi1, phi, gamma, p, log = TRUE)

  ProbOccNextTime <- psi1
  ll <- 0
  nyears <- dim(x)[1]
  if (nyears >= 1) {
    for (t in 1:nyears) {
      if (nrep[t] > 0) {
        numObs <- sum(x[t,1:nrep[t]])
        if (numObs < 0) {
          print("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
          stop("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
        }
        ProbOccAndCount <- ProbOccNextTime *
            exp(sum(dbinom(x[t,1:nrep[t]],
                           size = 1, p = p[t,1:nrep[t]], log = 1)))
        ProbUnoccAndCount <- (1-ProbOccNextTime) * (numObs == 0)
        ProbCount <- ProbOccAndCount + ProbUnoccAndCount
        ProbOccGivenCount <- ProbOccAndCount / ProbCount
        ll <- ll + log(ProbCount)
        if(t < nyears)
            ProbOccNextTime <- ProbOccGivenCount * phi[t] +
                (1-ProbOccGivenCount) * gamma
      } else {
        ## If there were no observations in year t,
        ## simply propagate probability of occupancy forward
        if (t < nyears)
            ProbOccNextTime <- ProbOccNextTime * phi[t] +
                (1-ProbOccNextTime) * gamma
      }
    }
  }

  lCorrectProbX <- ll

  expect_equal(exp(lCorrectProbX), probX)
  expect_equal(lCorrectProbX, lProbX)

  CdDynOcc <- compileNimble(dDynOcc_vs)
  CprobX <- CdDynOcc(x, nrep, psi1, phi, gamma, p, log = FALSE)
  expect_equal(CprobX, probX)

  ClProbX <- CdDynOcc(x, nrep, psi1, phi, gamma, p, log = TRUE)
  expect_equal(ClProbX, lProbX)



  nc <- nimbleCode({

    x[1:2, 1:5] ~ dDynOcc_vs(nrep[1:2], psi1, phi[1:2], gamma, p[1:2,1:5])

    psi1 ~ dunif(0,1)

    for (i in 1:2) {
      phi[i] ~ dunif(0,1)
    }
    gamma ~ dunif(0,1)

    for (i in 1:2) {
      for (j in 1:5) {
        p[i,j] ~ dunif(0,1)
      }
    }
  })

  m <- nimbleModel(nc, data = list(x = x, nrep = nrep),
                   # constants = list(nrep = nrep),
                   inits = list(p = p, phi = phi,
                                psi1 = psi1, gamma = gamma))

  m$calculate()
  MlProbX <- m$getLogProb("x")
  expect_equal(MlProbX, lProbX)

  cm <- compileNimble(m, showCompilerOutput = TRUE)
  cm$calculate()
  CMlProbX <- cm$getLogProb("x")
  expect_equal(CMlProbX, lProbX)

  set.seed(2468)
  cm$simulate('x')
  expect_equal(cm$x, matrix(c(0, 0, 1, 0, 1, 1, 1, 1, 1, 0), nrow = 2))
})

test_that("dDynOcc_sv works", {
  x <- matrix(c(0,0,1,0,1,
                1,1,1,1,0), nrow = 2)
  nrep <- c(5, 5)
  psi1 <- 0.7
  phi <- 0.8
  gamma <- c(0.1, 0.2)
  p <- matrix(rep(0.8, 10), nrow = 2)

  probX <- dDynOcc_sv(x, nrep, psi1, phi, gamma, p, log = FALSE)
  lProbX <- dDynOcc_sv(x, nrep, psi1, phi, gamma, p, log = TRUE)

  ProbOccNextTime <- psi1
  ll <- 0
  nyears <- dim(x)[1]
  if (nyears >= 1) {
    for (t in 1:nyears) {
      if (nrep[t] > 0) {
        numObs <- sum(x[t,1:nrep[t]])
        if (numObs < 0) {
          print("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
          stop("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
        }
        ProbOccAndCount <- ProbOccNextTime *
            exp(sum(dbinom(x[t,1:nrep[t]],
                           size = 1, p = p[t,1:nrep[t]], log = 1)))
        ProbUnoccAndCount <- (1-ProbOccNextTime) * (numObs == 0)
        ProbCount <- ProbOccAndCount + ProbUnoccAndCount
        ProbOccGivenCount <- ProbOccAndCount / ProbCount
        ll <- ll + log(ProbCount)
        if(t < nyears)
            ProbOccNextTime <- ProbOccGivenCount * phi +
                (1-ProbOccGivenCount) * gamma[t]
      } else {
        ## If there were no observations in year t,
        ## simply propagate probability of occupancy forward
        if (t < nyears)
            ProbOccNextTime <- ProbOccNextTime * phi +
                (1-ProbOccNextTime) * gamma[t]
      }
    }
  }

  lCorrectProbX <- ll

  expect_equal(exp(lCorrectProbX), probX)
  expect_equal(lCorrectProbX, lProbX)

  CdDynOcc <- compileNimble(dDynOcc_sv)
  CprobX <- CdDynOcc(x, nrep, psi1, phi, gamma, p, log = FALSE)
  expect_equal(CprobX, probX)

  ClProbX <- CdDynOcc(x, nrep, psi1, phi, gamma, p, log = TRUE)
  expect_equal(ClProbX, lProbX)



  nc <- nimbleCode({

    x[1:2, 1:5] ~ dDynOcc_sv(nrep[1:2], psi1, phi, gamma[1:2], p[1:2,1:5])

    psi1 ~ dunif(0,1)

    phi ~ dunif(0,1)
    for (i in 1:2) {
      gamma[i] ~ dunif(0,1)
    }


    for (i in 1:2) {
      for (j in 1:5) {
        p[i,j] ~ dunif(0,1)
      }
    }
  })

  m <- nimbleModel(nc, data = list(x = x, nrep = nrep),
                   # constants = list(nrep = nrep),
                   inits = list(p = p, phi = phi,
                                psi1 = psi1, gamma = gamma))

  m$calculate()
  MlProbX <- m$getLogProb("x")
  expect_equal(MlProbX, lProbX)

  cm <- compileNimble(m, showCompilerOutput = TRUE)
  cm$calculate()
  CMlProbX <- cm$getLogProb("x")
  expect_equal(CMlProbX, lProbX)

  set.seed(2468)
  cm$simulate('x')
  expect_equal(cm$x, matrix(c(0, 0, 1, 0, 1, 1, 1, 1, 1, 0), nrow = 2))
})


test_that("dDynOcc_ss works", {
  x <- matrix(c(0,0,1,0,1,
                1,1,1,1,0), nrow = 2)
  nrep <- c(5, 5)
  psi1 <- 0.7
  phi <- 0.8
  gamma <- 0.5
  p <- matrix(rep(0.8, 10), nrow = 2)

  probX <- dDynOcc_ss(x, nrep, psi1, phi, gamma, p, log = FALSE)
  lProbX <- dDynOcc_ss(x, nrep, psi1, phi, gamma, p, log = TRUE)

  ProbOccNextTime <- psi1
  ll <- 0
  nyears <- dim(x)[1]
  if (nyears >= 1) {
    for (t in 1:nyears) {
      if (nrep[t] > 0) {
        numObs <- sum(x[t,1:nrep[t]])
        if (numObs < 0) {
          print("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
          stop("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
        }
        ProbOccAndCount <- ProbOccNextTime *
            exp(sum(dbinom(x[t,1:nrep[t]],
                           size = 1, p = p[t,1:nrep[t]], log = 1)))
        ProbUnoccAndCount <- (1-ProbOccNextTime) * (numObs == 0)
        ProbCount <- ProbOccAndCount + ProbUnoccAndCount
        ProbOccGivenCount <- ProbOccAndCount / ProbCount
        ll <- ll + log(ProbCount)
        if(t < nyears)
            ProbOccNextTime <- ProbOccGivenCount * phi +
                (1-ProbOccGivenCount) * gamma
      } else {
        ## If there were no observations in year t,
        ## simply propagate probability of occupancy forward
        if (t < nyears)
            ProbOccNextTime <- ProbOccNextTime * phi +
                (1-ProbOccNextTime) * gamma
      }
    }
  }

  lCorrectProbX <- ll

  expect_equal(exp(lCorrectProbX), probX)
  expect_equal(lCorrectProbX, lProbX)

  CdDynOcc <- compileNimble(dDynOcc_ss)
  CprobX <- CdDynOcc(x, nrep, psi1, phi, gamma, p, log = FALSE)
  expect_equal(CprobX, probX)

  ClProbX <- CdDynOcc(x, nrep, psi1, phi, gamma, p, log = TRUE)
  expect_equal(ClProbX, lProbX)

  nc <- nimbleCode({

    x[1:2, 1:5] ~ dDynOcc_ss(nrep[1:2], psi1, phi, gamma, p[1:2,1:5])

    psi1 ~ dunif(0,1)

    phi ~ dunif(0,1)
    gamma ~ dunif(0,1)

    for (i in 1:2) {
      for (j in 1:5) {
        p[i,j] ~ dunif(0,1)
      }
    }
  })

  m <- nimbleModel(nc, data = list(x = x, nrep = nrep),
                   # constants = list(nrep = nrep),
                   inits = list(p = p, phi = phi,
                                psi1 = psi1, gamma = gamma))

  m$calculate()
  MlProbX <- m$getLogProb("x")
  expect_equal(MlProbX, lProbX)

  cm <- compileNimble(m, showCompilerOutput = TRUE)
  cm$calculate()
  CMlProbX <- cm$getLogProb("x")
  expect_equal(CMlProbX, lProbX)

  set.seed(2468)
  cm$simulate('x')
  expect_equal(cm$x, matrix(c(0, 0, 1, 0, 1, 1, 1, 1, 1, 0), nrow = 2))
})


test_that("Case errors in compiled dDynOcc_** work", {
  CdDynOcc_ss <- compileNimble(dDynOcc_ss)
  CdDynOcc_sv <- compileNimble(dDynOcc_sv)
  CdDynOcc_vs <- compileNimble(dDynOcc_vs)
  CdDynOcc_vv <- compileNimble(dDynOcc_vv)


  x <- matrix(c(0,0,1,0,1,
                1,1,1,1,0), nrow = 2)
  nrep <- c(5, 5)
  psi1 <- 0.7
  phi <- 0.8
  gamma <- 0.5
  p <- matrix(rep(0.8, 10), nrow = 2)

  expect_error(CdDynOcc_sv(x, nrep, psi1, phi, gamma, p, log = FALSE))
  expect_error(CdDynOcc_vs(x, nrep, psi1, phi, gamma, p, log = FALSE))
  expect_error(CdDynOcc_vv(x, nrep, psi1, phi, gamma, p, log = FALSE))

  phi <- c(0.8, 0.5)
  expect_error(CdDynOcc_sv(x, nrep, psi1, phi, gamma, p, log = FALSE))
  # expect_error(CdDynOcc_ss(x, nrep, psi1, phi, gamma, p, log = FALSE))
  expect_error(CdDynOcc_vv(x, nrep, psi1, phi, gamma, p, log = FALSE))

  gamma <- c(0.2, 0.1)
  # expect_error(CdDynOcc_sv(x, nrep, psi1, phi, gamma, p, log = FALSE))
  # expect_error(CdDynOcc_ss(x, nrep, psi1, phi, gamma, p, log = FALSE))
  # expect_error(CdDynOcc_vs(x, nrep, psi1, phi, gamma, p, log = FALSE))

  phi <- 0.3
  expect_error(CdDynOcc_vv(x, nrep, psi1, phi, gamma, p, log = FALSE))
  # expect_error(CdDynOcc_ss(x, nrep, psi1, phi, gamma, p, log = FALSE))
  # expect_error(CdDynOcc_vs(x, nrep, psi1, phi, gamma, p, log = FALSE))

  p <- p[1:2, 1:4]

  expect_error(CdDynOcc_sv(x, nrep, psi1, phi, gamma, p, log = FALSE))
  phi <- c(0.3, 0.5)
  expect_error(CdDynOcc_vv(x, nrep, psi1, phi, gamma, p, log = FALSE))
  gamma <- 0.5
  expect_error(CdDynOcc_vs(x, nrep, psi1, phi, gamma, p, log = FALSE))
  phi <- 0.3
  expect_error(CdDynOcc_ss(x, nrep, psi1, phi, gamma, p, log = FALSE))

})


