
context("Testing dDynOcc-related functions.")

test_that("dDynOcc_vvm works", {
  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)
  psi1 <- 0.7
  phi <- c(0.4, 0.4, 0.1)
  gamma <- c(0.4, 0.2, 0.1)
  p <- matrix(rep(0.8, 20), nrow = 4)

  probX <- dDynOcc_vvm(x, psi1, phi, gamma, p, start, end, log = FALSE)
  lProbX <- dDynOcc_vvm(x, psi1, phi, gamma, p, start, end, log = TRUE)

  ProbOccNextTime <- psi1
  ll <- 0
  nyears <- dim(x)[1]
  if (nyears >= 1) {
    for (t in 1:nyears) {
      if (end[t] - start[t] > 0) {
        numObs <- sum(x[t,start[t]:end[t]])
        if (numObs < 0) {
          print("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
          stop("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
        }
        ProbOccAndCount <- ProbOccNextTime *
            exp(sum(dbinom(x[t,start[t]:end[t]],
                           size = 1, prob = p[t,start[t]:end[t]], log = 1)))
        ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
        ProbCount <- ProbOccAndCount + ProbUnoccAndCount
        ProbOccGivenCount <- ProbOccAndCount / ProbCount
        ll <- ll + log(ProbCount)
        if (t < nyears)
            ProbOccNextTime <- ProbOccGivenCount * phi[t] +
                (1 - ProbOccGivenCount) * gamma[t]
      } else {
        ## If there were no observations in year t,
        ## simply propagate probability of occupancy forward
        if (t < nyears)
            ProbOccNextTime <- ProbOccNextTime *  phi[t] +
                (1 - ProbOccNextTime) * gamma[t]
      }
    }
  }

  lCorrectProbX <- ll

  expect_equal(exp(lCorrectProbX), probX)
  expect_equal(lCorrectProbX, lProbX)

  CdDynOcc_vvm <- compileNimble(dDynOcc_vvm)
  CprobX <- CdDynOcc_vvm(x, psi1, phi, gamma, p, start, end, log = FALSE)
  expect_equal(CprobX, probX)

  ClProbX <- CdDynOcc_vvm(x, psi1, phi, gamma, p, start, end, log = TRUE)
  expect_equal(ClProbX, lProbX)

  nc <- nimbleCode({

    x[1:4, 1:5] ~ dDynOcc_vvm(psi1, phi[1:3],
                              gamma[1:3], p[1:4,1:5],
                              start[1:4], end[1:4])

    psi1 ~ dunif(0,1)

    for (i in 1:3) {
      phi[i] ~ dunif(0,1)
      gamma[i] ~ dunif(0,1)
    }

    for (i in 1:4) {
      for (j in 1:5) {
        p[i,j] ~ dunif(0,1)
      }
    }
  })

  m <- nimbleModel(code = nc, data = list(x = x),
                   constants = list(start = start, end = end),
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


test_that("dDynOcc_vsm works", {
  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)
  psi1 <- 0.7
  phi <- c(0.4, 0.4, 0.1)
  gamma <- 0.5
  p <- matrix(rep(0.8, 20), nrow = 4)

  probX <- dDynOcc_vsm(x, psi1, phi, gamma, p, start, end, log = FALSE)
  lProbX <- dDynOcc_vsm(x, psi1, phi, gamma, p, start, end, log = TRUE)

  ProbOccNextTime <- psi1
  ll <- 0
  nyears <- dim(x)[1]
  if (nyears >= 1) {
    if (nyears >= 1) {
      for (t in 1:nyears) {
        if (end[t] - start[t] + 1 > 0) {
          numObs <- sum(x[t,start[t]:end[t]])
          if (numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,start[t]:end[t]],
                             size = 1, prob = p[t,start[t]:end[t]], log = 1)))
          ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
          ProbCount <- ProbOccAndCount + ProbUnoccAndCount
          ProbOccGivenCount <- ProbOccAndCount / ProbCount
          ll <- ll + log(ProbCount)
          if (t < nyears)
              ProbOccNextTime <- ProbOccGivenCount * phi[t] +
                  (1 - ProbOccGivenCount) * gamma
        } else {
          ## If there were no observations in year t,
          ## simply propagate probability of occupancy forward
          if (t < nyears)
              ProbOccNextTime <- ProbOccNextTime *  phi[t] +
                  (1 - ProbOccNextTime) * gamma
        }
      }
    }
  }

  lCorrectProbX <- ll

  expect_equal(exp(lCorrectProbX), probX)
  expect_equal(lCorrectProbX, lProbX)

  CdDynOcc_vsm <- compileNimble(dDynOcc_vsm)
  CprobX <- CdDynOcc_vsm(x, psi1, phi, gamma, p, start, end, log = FALSE)
  expect_equal(CprobX, probX)

  ClProbX <- CdDynOcc_vsm(x, psi1, phi, gamma, p, start, end, log = TRUE)
  expect_equal(ClProbX, lProbX)

  nc <- nimbleCode({

    x[1:4, 1:5] ~ dDynOcc_vsm(psi1, phi[1:3], gamma, p[1:4,1:5],
                              start[1:4], end[1:4])

    psi1 ~ dunif(0,1)

    for (i in 1:3) {
      phi[i] ~ dunif(0,1)
    }
    gamma ~ dunif(0,1)

    for (i in 1:4) {
      for (j in 1:5) {
        p[i,j] ~ dunif(0,1)
      }
    }
  })

  m <- nimbleModel(nc, data = list(x = x),
                   constants = list(start = start, end = end),
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

test_that("dDynOcc_svm works", {
  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)
  psi1 <- 0.7
  phi <- 0.5
  gamma <- c(0.4, 0.4, 0.1)
  p <- matrix(rep(0.8, 20), nrow = 4)

  probX <- dDynOcc_svm(x, psi1, phi, gamma, p, start, end, log = FALSE)
  lProbX <- dDynOcc_svm(x, psi1, phi, gamma, p, start, end, log = TRUE)

  ProbOccNextTime <- psi1
  ll <- 0
  nyears <- dim(x)[1]
    if(nyears >= 1) {
      for(t in 1:nyears) {
        if (end[t] - start[t] + 1 > 0) {
          numObs <- sum(x[t,start[t]:end[t]])
          if (numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,start[t]:end[t]],
                             size = 1, prob = p[t,start[t]:end[t]], log = 1)))
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

  CdDynOcc_svm <- compileNimble(dDynOcc_svm)
  CprobX <- CdDynOcc_svm(x, psi1, phi, gamma, p, start, end, log = FALSE)
  expect_equal(CprobX, probX)

  ClProbX <- CdDynOcc_svm(x, psi1, phi, gamma, p, start, end, log = TRUE)
  expect_equal(ClProbX, lProbX)



  nc <- nimbleCode({

    x[1:4, 1:5] ~ dDynOcc_svm(psi1, phi, gamma[1:3], p[1:4,1:5],
                             start[1:4], end[1:4])

    psi1 ~ dunif(0,1)

    phi ~ dunif(0,1)
    for (i in 1:3) {
      gamma[i] ~ dunif(0,1)
    }


    for (i in 1:4) {
      for (j in 1:5) {
        p[i,j] ~ dunif(0,1)
      }
    }
  })

  m <- nimbleModel(nc, data = list(x = x),
                   constants = list(start = start, end = end),
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


test_that("dDynOcc_ssm works", {
  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)
  psi1 <- 0.7
  phi <- 0.5
  gamma <- 0.2
  p <- matrix(rep(0.8, 20), nrow = 4)

  probX <- dDynOcc_ssm(x, psi1, phi, gamma, p, start, end, log = FALSE)
  lProbX <- dDynOcc_ssm(x, psi1, phi, gamma, p, start, end, log = TRUE)

  ProbOccNextTime <- psi1
  ll <- 0
  nyears <- dim(x)[1]
  if (nyears >= 1) {
    for (t in 1:nyears) {
      if (end[t] - start[t] + 1 > 0) {
        numObs <- sum(x[t,start[t]:end[t]])
        if (numObs < 0) {
          print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
        }
        ProbOccAndCount <- ProbOccNextTime *
            exp(sum(dbinom(x[t,start[t]:end[t]],
                           size = 1, prob = p[t,start[t]:end[t]], log = 1)))
        ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
        ProbCount <- ProbOccAndCount + ProbUnoccAndCount
        ProbOccGivenCount <- ProbOccAndCount / ProbCount
        ll <- ll + log(ProbCount)
        if (t < nyears)
            ProbOccNextTime <- ProbOccGivenCount * phi +
                (1 - ProbOccGivenCount) * gamma
      } else {
        ## If there were no observations in year t,
        ## simply propagate probability of occupancy forward
        if (t < nyears)
            ProbOccNextTime <- ProbOccNextTime * phi +
                (1 - ProbOccNextTime) * gamma
      }
    }
  }

  lCorrectProbX <- ll

  expect_equal(exp(lCorrectProbX), probX)
  expect_equal(lCorrectProbX, lProbX)

  CdDynOcc_ssm <- compileNimble(dDynOcc_ssm)
  CprobX <- CdDynOcc_ssm(x, psi1, phi, gamma, p, start, end, log = FALSE)
  expect_equal(CprobX, probX)

  ClProbX <- CdDynOcc_ssm(x, psi1, phi, gamma, p, start, end, log = TRUE)
  expect_equal(ClProbX, lProbX)

  nc <- nimbleCode({

    x[1:4, 1:5] ~ dDynOcc_ssm(psi1, phi, gamma, p[1:4, 1:5],
                              start[1:4], end[1:4])

    psi1 ~ dunif(0,1)

    phi ~ dunif(0,1)
    gamma ~ dunif(0,1)

    for (i in 1:2) {
      for (j in 1:5) {
        p[i,j] ~ dunif(0,1)
      }
    }
  })

  m <- nimbleModel(code = nc, data = list(x = x),
                   constants = list(start = start, end = end),
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


test_that("Case errors in compiled dDynOcc_**m work", {
  CdDynOcc_ssm <- compileNimble(dDynOcc_ssm)
  CdDynOcc_svm <- compileNimble(dDynOcc_svm)
  CdDynOcc_vsm <- compileNimble(dDynOcc_vsm)
  CdDynOcc_vvm <- compileNimble(dDynOcc_vvm)


  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  nrep <- c(5, 5)
  psi1 <- 0.7
  phi <- 0.8
  gamma <- 0.5
  p <- matrix(rep(0.8, 20), nrow = 4)

  expect_error(CdDynOcc_svm(x, psi1, phi, gamma, p, start, end, log = FALSE))
  expect_error(CdDynOcc_vsm(x, psi1, phi, gamma, p, start, end, log = FALSE))
  expect_error(CdDynOcc_vvm(x, psi1, phi, gamma, p, start, end, log = FALSE))

  phi <- c(0.8, 0.5, 0.2)
  expect_error(CdDynOcc_svm(x, psi1, phi, gamma, p, start, end, log = FALSE))
  expect_error(CdDynOcc_vvm(x, psi1, phi, gamma, p, start, end, log = FALSE))

  gamma <- c(0.2, 0.5, 0.5)

  phi <- 0.3
  expect_error(CdDynOcc_vvm(x, psi1, phi, gamma, p, start, end, log = FALSE))

  p <- p[1:2, 1:4]

  expect_error(CdDynOcc_svm(x, psi1, phi, gamma, p, start, end, log = FALSE))
  phi <- c(0.3, 0.5, 0.3)
  expect_error(CdDynOcc_vvm(x, psi1, phi, gamma, p, start, end, log = FALSE))
  gamma <- 0.5
  expect_error(CdDynOcc_vsm(x, psi1, phi, gamma, p, start, end, log = FALSE))
  phi <- 0.3
  expect_error(CdDynOcc_ssm(x, psi1, phi, gamma, p, start, end, log = FALSE))

})




test_that("dDynOcc_vvv works", {
  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)
  psi1 <- 0.7
  phi <- c(0.4, 0.4, 0.1)
  gamma <- c(0.4, 0.2, 0.1)
  p <- rep(0.8, 4)

  probX <- dDynOcc_vvv(x, psi1, phi, gamma, p, start, end, log = FALSE)
  lProbX <- dDynOcc_vvv(x, psi1, phi, gamma, p, start, end, log = TRUE)

  ProbOccNextTime <- psi1
  ll <- 0
  nyears <- dim(x)[1]
  if (nyears >= 1) {
    for (t in 1:nyears) {
      if (end[t] - start[t] > 0) {
        numObs <- sum(x[t,start[t]:end[t]])
        if (numObs < 0) {
          print("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
          stop("Error in dDynamicOccupancy: numObs < 0 but nrep[t] > 0\n")
        }
        ProbOccAndCount <- ProbOccNextTime *
            exp(sum(dbinom(x[t,start[t]:end[t]],
                           size = 1, prob = p[t], log = 1)))
        ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
        ProbCount <- ProbOccAndCount + ProbUnoccAndCount
        ProbOccGivenCount <- ProbOccAndCount / ProbCount
        ll <- ll + log(ProbCount)
        if (t < nyears)
            ProbOccNextTime <- ProbOccGivenCount * phi[t] +
                (1 - ProbOccGivenCount) * gamma[t]
      } else {
        ## If there were no observations in year t,
        ## simply propagate probability of occupancy forward
        if (t < nyears)
            ProbOccNextTime <- ProbOccNextTime *  phi[t] +
                (1 - ProbOccNextTime) * gamma[t]
      }
    }
  }

  lCorrectProbX <- ll

  expect_equal(exp(lCorrectProbX), probX)
  expect_equal(lCorrectProbX, lProbX)

  CdDynOcc_vvv <- compileNimble(dDynOcc_vvv)
  CprobX <- CdDynOcc_vvv(x, psi1, phi, gamma, p, start, end, log = FALSE)
  expect_equal(CprobX, probX)

  ClProbX <- CdDynOcc_vvv(x, psi1, phi, gamma, p, start, end, log = TRUE)
  expect_equal(ClProbX, lProbX)

  nc <- nimbleCode({

    x[1:4, 1:5] ~ dDynOcc_vvv(psi1 = psi1, phi = phi[1:3],
                              gamma = gamma[1:3], p = p[1:4],
                              start = start[1:4], end = end[1:4])

    psi1 ~ dunif(0,1)

    for (i in 1:3) {
      phi[i] ~ dunif(0,1)
      gamma[i] ~ dunif(0,1)
    }

    for (i in 1:4) {
      p[i] ~ dunif(0,1)
    }
  })

  m <- nimbleModel(code = nc, data = list(x = x),
                   constants = list(start = start, end = end),
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


test_that("dDynOcc_vsv works", {
  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)
  psi1 <- 0.7
  phi <- c(0.4, 0.4, 0.1)
  gamma <- 0.5
  p <- rep(0.8, 4)

  probX <- dDynOcc_vsv(x, psi1, phi, gamma, p, start, end, log = FALSE)
  lProbX <- dDynOcc_vsv(x, psi1, phi, gamma, p, start, end, log = TRUE)

  ProbOccNextTime <- psi1
  ll <- 0
  nyears <- dim(x)[1]
  if (nyears >= 1) {
    if (nyears >= 1) {
      for (t in 1:nyears) {
        if (end[t] - start[t] + 1 > 0) {
          numObs <- sum(x[t,start[t]:end[t]])
          if (numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,start[t]:end[t]],
                             size = 1, prob = p[t], log = 1)))
          ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
          ProbCount <- ProbOccAndCount + ProbUnoccAndCount
          ProbOccGivenCount <- ProbOccAndCount / ProbCount
          ll <- ll + log(ProbCount)
          if (t < nyears)
              ProbOccNextTime <- ProbOccGivenCount * phi[t] +
                  (1 - ProbOccGivenCount) * gamma
        } else {
          ## If there were no observations in year t,
          ## simply propagate probability of occupancy forward
          if (t < nyears)
              ProbOccNextTime <- ProbOccNextTime *  phi[t] +
                  (1 - ProbOccNextTime) * gamma
        }
      }
    }
  }

  lCorrectProbX <- ll

  expect_equal(exp(lCorrectProbX), probX)
  expect_equal(lCorrectProbX, lProbX)

  CdDynOcc_vsv <- compileNimble(dDynOcc_vsv)
  CprobX <- CdDynOcc_vsv(x, psi1, phi, gamma, p, start, end, log = FALSE)
  expect_equal(CprobX, probX)

  ClProbX <- CdDynOcc_vsv(x, psi1, phi, gamma, p, start, end, log = TRUE)
  expect_equal(ClProbX, lProbX)

  nc <- nimbleCode({

    x[1:4, 1:5] ~ dDynOcc_vsv(psi1, phi[1:3], gamma, p[1:4],
                              start[1:4], end[1:4])

    psi1 ~ dunif(0,1)

    for (i in 1:3) {
      phi[i] ~ dunif(0,1)
    }
    gamma ~ dunif(0,1)

    for (i in 1:4) {
      p[i] ~ dunif(0,1)
    }
  })

  m <- nimbleModel(nc, data = list(x = x),
                   constants = list(start = start, end = end),
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

test_that("dDynOcc_svv works", {
  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)
  psi1 <- 0.7
  phi <- 0.5
  gamma <- c(0.4, 0.4, 0.1)
  p <- rep(0.8, 4)

  probX <- dDynOcc_svv(x, psi1, phi, gamma, p, start, end, log = FALSE)
  lProbX <- dDynOcc_svv(x, psi1, phi, gamma, p, start, end, log = TRUE)

  ProbOccNextTime <- psi1
  ll <- 0
  nyears <- dim(x)[1]
    if(nyears >= 1) {
      for(t in 1:nyears) {
        if (end[t] - start[t] + 1 > 0) {
          numObs <- sum(x[t,start[t]:end[t]])
          if (numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,start[t]:end[t]],
                             size = 1, prob = p[t], log = 1)))
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

  CdDynOcc_svv <- compileNimble(dDynOcc_svv)
  CprobX <- CdDynOcc_svv(x, psi1, phi, gamma, p, start, end, log = FALSE)
  expect_equal(CprobX, probX)

  ClProbX <- CdDynOcc_svv(x, psi1, phi, gamma, p, start, end, log = TRUE)
  expect_equal(ClProbX, lProbX)



  nc <- nimbleCode({

    x[1:4, 1:5] ~ dDynOcc_svv(psi1, phi, gamma[1:3], p[1:4],
                             start[1:4], end[1:4])

    psi1 ~ dunif(0,1)

    phi ~ dunif(0,1)
    for (i in 1:3) {
      gamma[i] ~ dunif(0,1)
    }


    for (i in 1:4) {
      p[i,j] ~ dunif(0,1)
    }
  })

  m <- nimbleModel(nc, data = list(x = x),
                   constants = list(start = start, end = end),
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


test_that("dDynOcc_ssv works", {
  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)
  psi1 <- 0.7
  phi <- 0.5
  gamma <- 0.2
  p <- rep(0.8, 4)

  probX <- dDynOcc_ssv(x, psi1, phi, gamma, p, start, end, log = FALSE)
  lProbX <- dDynOcc_ssv(x, psi1, phi, gamma, p, start, end, log = TRUE)

  ProbOccNextTime <- psi1
  ll <- 0
  nyears <- dim(x)[1]
  if (nyears >= 1) {
    for (t in 1:nyears) {
      if (end[t] - start[t] + 1 > 0) {
        numObs <- sum(x[t,start[t]:end[t]])
        if (numObs < 0) {
          print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
        }
        ProbOccAndCount <- ProbOccNextTime *
            exp(sum(dbinom(x[t,start[t]:end[t]],
                           size = 1, prob = p[t], log = 1)))
        ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
        ProbCount <- ProbOccAndCount + ProbUnoccAndCount
        ProbOccGivenCount <- ProbOccAndCount / ProbCount
        ll <- ll + log(ProbCount)
        if (t < nyears)
            ProbOccNextTime <- ProbOccGivenCount * phi +
                (1 - ProbOccGivenCount) * gamma
      } else {
        ## If there were no observations in year t,
        ## simply propagate probability of occupancy forward
        if (t < nyears)
            ProbOccNextTime <- ProbOccNextTime * phi +
                (1 - ProbOccNextTime) * gamma
      }
    }
  }

  lCorrectProbX <- ll

  expect_equal(exp(lCorrectProbX), probX)
  expect_equal(lCorrectProbX, lProbX)

  CdDynOcc_ssv <- compileNimble(dDynOcc_ssv)
  CprobX <- CdDynOcc_ssv(x, psi1, phi, gamma, p, start, end, log = FALSE)
  expect_equal(CprobX, probX)

  ClProbX <- CdDynOcc_ssv(x, psi1, phi, gamma, p, start, end, log = TRUE)
  expect_equal(ClProbX, lProbX)

  nc <- nimbleCode({

    x[1:4, 1:5] ~ dDynOcc_ssv(psi1, phi, gamma, p[1:4],
                              start[1:4], end[1:4])

    psi1 ~ dunif(0,1)

    phi ~ dunif(0,1)
    gamma ~ dunif(0,1)

    for (i in 1:4) {
      p[i] ~ dunif(0,1)
    }
  })

  m <- nimbleModel(code = nc, data = list(x = x),
                   constants = list(start = start, end = end),
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


test_that("Case errors in compiled dDynOcc_**v work", {
  CdDynOcc_ssv <- compileNimble(dDynOcc_ssv)
  CdDynOcc_svv <- compileNimble(dDynOcc_svv)
  CdDynOcc_vsv <- compileNimble(dDynOcc_vsv)
  CdDynOcc_vvv <- compileNimble(dDynOcc_vvv)


  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  nrep <- c(5, 5)
  psi1 <- 0.7
  phi <- 0.8
  gamma <- 0.5
  p <- rep(0.8, 4)

  expect_error(CdDynOcc_svv(x, psi1, phi, gamma, p, start, end, log = FALSE))
  expect_error(CdDynOcc_vsv(x, psi1, phi, gamma, p, start, end, log = FALSE))
  expect_error(CdDynOcc_vvv(x, psi1, phi, gamma, p, start, end, log = FALSE))

  phi <- c(0.8, 0.5, 0.2)
  expect_error(CdDynOcc_svv(x, psi1, phi, gamma, p, start, end, log = FALSE))
  expect_error(CdDynOcc_vvv(x, psi1, phi, gamma, p, start, end, log = FALSE))

  gamma <- c(0.2, 0.5, 0.5)

  phi <- 0.3
  expect_error(CdDynOcc_vvv(x, psi1, phi, gamma, p, start, end, log = FALSE))

  p <- p[1:2]

  expect_error(CdDynOcc_svv(x, psi1, phi, gamma, p, start, end, log = FALSE))
  phi <- c(0.3, 0.5, 0.3)
  expect_error(CdDynOcc_vvv(x, psi1, phi, gamma, p, start, end, log = FALSE))
  gamma <- 0.5
  expect_error(CdDynOcc_vsv(x, psi1, phi, gamma, p, start, end, log = FALSE))
  phi <- 0.3
  expect_error(CdDynOcc_ssv(x, psi1, phi, gamma, p, start, end, log = FALSE))

})


context("Testing dDynOcc-related functions.")

test_that("dDynOcc_vvs works", {
  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)
  psi1 <- 0.7
  phi <- c(0.4, 0.4, 0.1)
  gamma <- c(0.4, 0.2, 0.1)
  p <- 0.6

  probX <- dDynOcc_vvs(x, psi1, phi, gamma, p, start, end, log = FALSE)
  lProbX <- dDynOcc_vvs(x, psi1, phi, gamma, p, start, end, log = TRUE)

  ProbOccNextTime <- psi1
  ll <- 0
  nyears <- dim(x)[1]
  if (nyears >= 1) {
    for (t in 1:nyears) {
      if (end[t] - start[t] > 0) {
        numObs <- sum(x[t,start[t]:end[t]])
        if (numObs < 0) {
          print("Error in dDynamicOccupancy: numObs < 0 but nrep > 0\n")
          stop("Error in dDynamicOccupancy: numObs < 0 but nrep > 0\n")
        }
        ProbOccAndCount <- ProbOccNextTime *
            exp(sum(dbinom(x[t,start[t]:end[t]],
                           size = 1, prob = p, log = 1)))
        ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
        ProbCount <- ProbOccAndCount + ProbUnoccAndCount
        ProbOccGivenCount <- ProbOccAndCount / ProbCount
        ll <- ll + log(ProbCount)
        if (t < nyears)
            ProbOccNextTime <- ProbOccGivenCount * phi[t] +
                (1 - ProbOccGivenCount) * gamma[t]
      } else {
        ## If there were no observations in year t,
        ## simply propagate probability of occupancy forward
        if (t < nyears)
            ProbOccNextTime <- ProbOccNextTime *  phi[t] +
                (1 - ProbOccNextTime) * gamma[t]
      }
    }
  }

  lCorrectProbX <- ll

  expect_equal(exp(lCorrectProbX), probX)
  expect_equal(lCorrectProbX, lProbX)

  CdDynOcc_vvs <- compileNimble(dDynOcc_vvs)
  CprobX <- CdDynOcc_vvs(x, psi1, phi, gamma, p, start, end, log = FALSE)
  expect_equal(CprobX, probX)

  ClProbX <- CdDynOcc_vvs(x, psi1, phi, gamma, p, start, end, log = TRUE)
  expect_equal(ClProbX, lProbX)

  nc <- nimbleCode({

    x[1:4, 1:5] ~ dDynOcc_vvs(psi1 = psi1, phi = phi[1:3],
                              gamma = gamma[1:3], p = p,
                              start = start[1:4], end = end[1:4])

    psi1 ~ dunif(0,1)

    for (i in 1:3) {
      phi[i] ~ dunif(0,1)
      gamma[i] ~ dunif(0,1)
    }

    p ~ dunif(0,1)
  })

  m <- nimbleModel(code = nc, data = list(x = x),
                   constants = list(start = start, end = end),
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


test_that("dDynOcc_vss works", {
  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)
  psi1 <- 0.7
  phi <- c(0.4, 0.4, 0.1)
  gamma <- 0.5
  p <- 0.6

  probX <- dDynOcc_vss(x, psi1, phi, gamma, p, start, end, log = FALSE)
  lProbX <- dDynOcc_vss(x, psi1, phi, gamma, p, start, end, log = TRUE)

  ProbOccNextTime <- psi1
  ll <- 0
  nyears <- dim(x)[1]
  if (nyears >= 1) {
    if (nyears >= 1) {
      for (t in 1:nyears) {
        if (end[t] - start[t] + 1 > 0) {
          numObs <- sum(x[t,start[t]:end[t]])
          if (numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,start[t]:end[t]],
                             size = 1, prob = p, log = 1)))
          ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
          ProbCount <- ProbOccAndCount + ProbUnoccAndCount
          ProbOccGivenCount <- ProbOccAndCount / ProbCount
          ll <- ll + log(ProbCount)
          if (t < nyears)
              ProbOccNextTime <- ProbOccGivenCount * phi[t] +
                  (1 - ProbOccGivenCount) * gamma
        } else {
          ## If there were no observations in year t,
          ## simply propagate probability of occupancy forward
          if (t < nyears)
              ProbOccNextTime <- ProbOccNextTime *  phi[t] +
                  (1 - ProbOccNextTime) * gamma
        }
      }
    }
  }

  lCorrectProbX <- ll

  expect_equal(exp(lCorrectProbX), probX)
  expect_equal(lCorrectProbX, lProbX)

  CdDynOcc_vss <- compileNimble(dDynOcc_vss)
  CprobX <- CdDynOcc_vss(x, psi1, phi, gamma, p, start, end, log = FALSE)
  expect_equal(CprobX, probX)

  ClProbX <- CdDynOcc_vss(x, psi1, phi, gamma, p, start, end, log = TRUE)
  expect_equal(ClProbX, lProbX)

  nc <- nimbleCode({

    x[1:4, 1:5] ~ dDynOcc_vss(psi1, phi[1:3], gamma, p,
                              start[1:4], end[1:4])

    psi1 ~ dunif(0,1)

    for (i in 1:3) {
      phi[i] ~ dunif(0,1)
    }
    gamma ~ dunif(0,1)

    p ~ dunif(0,1)
  })

  m <- nimbleModel(nc, data = list(x = x),
                   constants = list(start = start, end = end),
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

test_that("dDynOcc_svs works", {
  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)
  psi1 <- 0.7
  phi <- 0.5
  gamma <- c(0.4, 0.4, 0.1)
  p <- 0.6

  probX <- dDynOcc_svs(x, psi1, phi, gamma, p, start, end, log = FALSE)
  lProbX <- dDynOcc_svs(x, psi1, phi, gamma, p, start, end, log = TRUE)

  ProbOccNextTime <- psi1
  ll <- 0
  nyears <- dim(x)[1]
    if(nyears >= 1) {
      for(t in 1:nyears) {
        if (end[t] - start[t] + 1 > 0) {
          numObs <- sum(x[t,start[t]:end[t]])
          if (numObs < 0) {
            print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
            stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          }
          ProbOccAndCount <- ProbOccNextTime *
              exp(sum(dbinom(x[t,start[t]:end[t]],
                             size = 1, prob = p, log = 1)))
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

  CdDynOcc_svs <- compileNimble(dDynOcc_svs)
  CprobX <- CdDynOcc_svs(x, psi1, phi, gamma, p, start, end, log = FALSE)
  expect_equal(CprobX, probX)

  ClProbX <- CdDynOcc_svs(x, psi1, phi, gamma, p, start, end, log = TRUE)
  expect_equal(ClProbX, lProbX)



  nc <- nimbleCode({

    x[1:4, 1:5] ~ dDynOcc_svs(psi1, phi, gamma[1:3], p,
                             start[1:4], end[1:4])

    psi1 ~ dunif(0,1)

    phi ~ dunif(0,1)
    for (i in 1:3) {
      gamma[i] ~ dunif(0,1)
    }

    p ~ dunif(0,1)
  })

  m <- nimbleModel(nc, data = list(x = x),
                   constants = list(start = start, end = end),
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


test_that("dDynOcc_sss works", {
  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)
  psi1 <- 0.7
  phi <- 0.5
  gamma <- 0.2
  p <- 0.8

  probX <- dDynOcc_sss(x, psi1, phi, gamma, p, start, end, log = FALSE)
  lProbX <- dDynOcc_sss(x, psi1, phi, gamma, p, start, end, log = TRUE)

  ProbOccNextTime <- psi1
  ll <- 0
  nyears <- dim(x)[1]
  if (nyears >= 1) {
    for (t in 1:nyears) {
      if (end[t] - start[t] + 1 > 0) {
        numObs <- sum(x[t,start[t]:end[t]])
        if (numObs < 0) {
          print("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
          stop("Error in dDynamicOccupancy: numObs < 0 but number of obs in start/end > 0\n")
        }
        ProbOccAndCount <- ProbOccNextTime *
            exp(sum(dbinom(x[t,start[t]:end[t]],
                           size = 1, prob = p, log = 1)))
        ProbUnoccAndCount <- (1 - ProbOccNextTime) * (numObs == 0)
        ProbCount <- ProbOccAndCount + ProbUnoccAndCount
        ProbOccGivenCount <- ProbOccAndCount / ProbCount
        ll <- ll + log(ProbCount)
        if (t < nyears)
            ProbOccNextTime <- ProbOccGivenCount * phi +
                (1 - ProbOccGivenCount) * gamma
      } else {
        ## If there were no observations in year t,
        ## simply propagate probability of occupancy forward
        if (t < nyears)
            ProbOccNextTime <- ProbOccNextTime * phi +
                (1 - ProbOccNextTime) * gamma
      }
    }
  }

  lCorrectProbX <- ll

  expect_equal(exp(lCorrectProbX), probX)
  expect_equal(lCorrectProbX, lProbX)

  CdDynOcc_sss <- compileNimble(dDynOcc_sss)
  CprobX <- CdDynOcc_sss(x, psi1, phi, gamma, p, start, end, log = FALSE)
  expect_equal(CprobX, probX)

  ClProbX <- CdDynOcc_sss(x, psi1, phi, gamma, p, start, end, log = TRUE)
  expect_equal(ClProbX, lProbX)

  nc <- nimbleCode({

    x[1:4, 1:5] ~ dDynOcc_sss(psi1, phi, gamma, p,
                              start[1:4], end[1:4])

    psi1 ~ dunif(0,1)

    phi ~ dunif(0,1)
    gamma ~ dunif(0,1)

    p ~ dunif(0,1)
  })

  m <- nimbleModel(code = nc, data = list(x = x),
                   constants = list(start = start, end = end),
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
  CdDynOcc_sss <- compileNimble(dDynOcc_sss)
  CdDynOcc_svs <- compileNimble(dDynOcc_svs)
  CdDynOcc_vss <- compileNimble(dDynOcc_vss)
  CdDynOcc_vvs <- compileNimble(dDynOcc_vvs)


  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  nrep <- c(5, 5)
  psi1 <- 0.7
  phi <- 0.8
  gamma <- 0.5
  p <- 0.8

  expect_error(CdDynOcc_svs(x, psi1, phi, gamma, p, start, end, log = FALSE))
  expect_error(CdDynOcc_vss(x, psi1, phi, gamma, p, start, end, log = FALSE))
  expect_error(CdDynOcc_vvs(x, psi1, phi, gamma, p, start, end, log = FALSE))

  phi <- c(0.8, 0.5, 0.2)
  expect_error(CdDynOcc_svs(x, psi1, phi, gamma, p, start, end, log = FALSE))
  expect_error(CdDynOcc_vvs(x, psi1, phi, gamma, p, start, end, log = FALSE))

  gamma <- c(0.2, 0.5, 0.5)

  phi <- 0.3
  expect_error(CdDynOcc_vvs(x, psi1, phi, gamma, p, start, end, log = FALSE))

})
