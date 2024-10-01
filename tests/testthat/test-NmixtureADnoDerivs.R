# Test the N-mixture distribution nimbleFunction.
library(nimbleEcology)
# -----------------------------------------------------------------------------
# 0. Load
# -----------------------------------------------------------------------------
#### 1. Test dNmixtureAD_v ####
test_that("dNmixtureAD_v works uncompiled",
{
  # Uncompiled calculation
  x <- c(1, 0, 1, 3, 0)
  lambda <- 8
  prob <- c(0.5, 0.3, 0.5, 0.4, 0.1)
  Nmin <- 0
  Nmax <- 250
  len <- 5

  probX <- dNmixtureAD_v(x, lambda, prob, Nmin, Nmax, len)
  # Manually calculate the correct answer
  correctProbX <- 0
  for (N in Nmin:Nmax) {
    correctProbX <- correctProbX + dpois(N, lambda) * prod(dbinom(x, N, prob))
  }
  expect_equal(probX, correctProbX)

  # Uncompiled log probability
  lProbX <- dNmixtureAD_v(x, lambda, prob, Nmin, Nmax, len, log = TRUE)
  lCorrectProbX <- log(correctProbX)
  expect_equal(lProbX, lCorrectProbX)

  # Other Nmin / Nmax
  Nmin <- 3
  for(Nmax in 3:6)  {
    dynProbX <- dNmixtureAD_v(x, lambda, prob, Nmin = Nmin, Nmax = Nmax, len)
    dynCorrectProbX <- 0
    for (N in Nmin:Nmax) {
      dynCorrectProbX <- dynCorrectProbX + dpois(N, lambda) * prod(dbinom(x, N, prob))
    }
    expect_equal(dynProbX, dynCorrectProbX)
  }
  Nmin <- 0
  Nmax <- 250
  # Compilation and compiled calculations
  call_dNmixtureAD_v <- nimbleFunction(
    name = "t1",
    run=function(x=double(1),
                 lambda=double(),
                 prob=double(1),
                 Nmin = double(0),
                 Nmax = double(0),
                 len=double(),
                 log = integer(0, default=0)) {
      return(dNmixtureAD_v(x,lambda,prob,Nmin,Nmax,len,log))
      returnType(double())
    }
  )
  CdNmixtureAD_v <- compileNimble(call_dNmixtureAD_v)
  CprobX <- CdNmixtureAD_v(x, lambda, prob, Nmin, Nmax, len)
  probX <- dNmixtureAD_v(x, lambda, prob, Nmin, Nmax, len)
  expect_equal(CprobX, probX)

  ClProbX <- CdNmixtureAD_v(x, lambda, prob, Nmin, Nmax, len, log = TRUE)
  lProbX <- dNmixtureAD_v(x, lambda, prob, Nmin, Nmax, len, log = TRUE)
  expect_equal(ClProbX, lProbX)

  expect_error(CdynProbX <- CdNmixtureAD_v(x, lambda, prob, Nmin = -1, Nmax = -1, len))

  # Use in Nimble model
  nc <- nimbleCode({
    x[1:5] ~ dNmixtureAD_v(lambda = lambda, prob = prob[1:5],
                         Nmin = Nmin, Nmax = Nmax,
                         len = len)
  })

  m <- nimbleModel(code = nc,
                   data = list(x = x),
                   inits = list(lambda = lambda,
                                prob = prob),
                   constants = list(Nmin = Nmin, Nmax = Nmax,
                                    len = len))
  m$calculate()
  MlProbX <- m$getLogProb("x")
  expect_equal(MlProbX, lProbX)

  # Compiled model
  cm <- compileNimble(m)
  cm$calculate()
  CMlProbX <- cm$getLogProb("x")
  expect_equal(CMlProbX, lProbX)

  # Missing values
  xna <- c(1, 0, NA, 3, 0)
  probXna <- dNmixtureAD_v(xna, lambda, prob, Nmin, Nmax, len)
  # Manually calculate the correct answer
  correctProbXna <- 0
  for (N in Nmin:Nmax) {
    correctProbXna <- correctProbXna + dpois(N, lambda) * prod(dbinom(xna, N, prob), na.rm=TRUE)
  }
  expect_equal(probXna, correctProbXna)

  CprobXna <- CdNmixtureAD_v(xna, lambda, prob, Nmin, Nmax, len)
  expect_equal(CprobXna, probXna)

  xna <- c(1,NA,NA,NA,NA)
  correctProbXna <- 0
  for (N in Nmin:Nmax) {
    correctProbXna <- correctProbXna + dpois(N, lambda) * prod(dbinom(xna, N, prob), na.rm=TRUE)
  }
  probXna <- CdNmixtureAD_v(xna, lambda, prob, Nmin, Nmax, len)
  expect_equal(probXna, correctProbXna)

  xna <- as.numeric(rep(NA, 5))
  expect_equal(CdNmixtureAD_v(xna, lambda, prob, Nmin, Nmax, len), 1)
  expect_equal(CdNmixtureAD_v(xna, lambda, prob, Nmin, Nmax, len, log=TRUE), 0)

  # Test imputing value for all NAs
  xNA <- c(NA, NA, NA, NA, NA)
  mNA <- nimbleModel(nc, data = list(x = xNA),
                     inits = list(lambda = lambda,
                                  prob = prob),
                     constants = list(Nmin = Nmin, Nmax = Nmax,
                                      len = len))


  mNAConf <- configureMCMC(mNA)
  mNAConf$addMonitors('x')
  mNA_MCMC <- buildMCMC(mNAConf)
  cmNA <- compileNimble(mNA, mNA_MCMC)

  set.seed(0)
  cmNA$mNA_MCMC$run(10)

  # Did the imputed values come back?
  expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x[2]"])))

  # Test simulation code
  nSim <- 10
  xSim <- array(NA, dim = c(nSim, len))
  set.seed(1)
  for (i in 1:nSim) {
    xSim[i,] <- rNmixtureAD_v(1, lambda, prob, Nmin, Nmax, len)
  }

  CrNmixtureAD_v <- compileNimble(rNmixtureAD_v)
  CxSim <- array(NA, dim = c(nSim, len))
  set.seed(1)
  for (i in 1:nSim) {
    CxSim[i,] <- CrNmixtureAD_v(1, lambda, prob, Nmin, Nmax, len)
  }
  expect_identical(xSim, CxSim)

  simNodes <- m$getDependencies(c('prob', 'lambda'), self = FALSE)
  mxSim <- array(NA, dim = c(nSim, len))
  set.seed(1)
  for(i in 1:nSim) {
    m$simulate(simNodes, includeData = TRUE)
    mxSim[i,] <- m$x
  }
  expect_identical(mxSim, xSim)

  CmxSim <- array(NA, dim = c(nSim, len))
  set.seed(1)
  for(i in 1:nSim) {
    cm$simulate(simNodes, includeData = TRUE)
    CmxSim[i,] <- cm$x
  }
  expect_identical(CmxSim, mxSim)
})

# -----------------------------------------------------------------------------
#### 2. Test dNmixtureAD_s ####
test_that("dNmixtureAD_s works",
          {
      # Uncompiled calculation
          x <- c(1, 0, 1, 3, 2)
          lambda <- 8
          prob <- 0.4
          Nmin <- 0
          Nmax <- 250
          len <- 5

          probX <- dNmixtureAD_s(x, lambda, prob, Nmin, Nmax, len)
      # Manually calculate the correct answer
          correctProbX <- 0
          for (N in Nmin:Nmax) {
            correctProbX <- correctProbX + dpois(N, lambda) * prod(dbinom(x, N, prob))
          }

          expect_equal(probX, correctProbX)

      # Uncompiled log probability
          lProbX <- dNmixtureAD_s(x, lambda, prob, Nmin, Nmax, len, log = TRUE)
          lCorrectProbX <- log(correctProbX)
          expect_equal(lProbX, lCorrectProbX)

      # Dynamic Nmin / Nmax
            expect_error(dynProbX <-
                           dNmixtureAD_s(x, lambda, prob, Nmin = -1, Nmax = -1, len))

  # Other Nmin / Nmax
  Nmin <- 3
  for(Nmax in 3:6)  {
    dynProbX <- dNmixtureAD_s(x, lambda, prob, Nmin = Nmin, Nmax = Nmax, len)
    dynCorrectProbX <- 0
    for (N in Nmin:Nmax) {
      dynCorrectProbX <- dynCorrectProbX + dpois(N, lambda) * prod(dbinom(x, N, prob))
    }
    expect_equal(dynProbX, dynCorrectProbX)
  }
  Nmin <- 0
  Nmax <- 250
  # Compilation and compiled calculations
  call_dNmixtureAD_s <- nimbleFunction(
    name = "t2",
    run=function(x=double(1),
                 lambda=double(),
                 prob=double(),
                 Nmin = double(0),
                 Nmax = double(0),
                 len=double(),
                 log = integer(0, default=0)) {
      return(dNmixtureAD_s(x,lambda,prob,Nmin,Nmax,len,log))
      returnType(double())
    }
  )


      # Compilation and compiled calculations
          CdNmixtureAD_s <- compileNimble(call_dNmixtureAD_s)
          CprobX <- CdNmixtureAD_s(x, lambda, prob, Nmin, Nmax, len)
          expect_equal(CprobX, probX)

          ClProbX <- CdNmixtureAD_s(x, lambda, prob, Nmin, Nmax, len, log = TRUE)
          expect_equal(ClProbX, lProbX)

          expect_error(CdynProbX <- CdNmixtureAD_s(x, lambda, prob, Nmin = -1, Nmax = -1, len))

      # Use in Nimble model
          nc <- nimbleCode({
            x[1:5] ~ dNmixtureAD_s(lambda = lambda, prob = prob,
                            Nmin = Nmin, Nmax = Nmax, len = len)
          })

          m <- nimbleModel(code = nc,
                           data = list(x = x),
                           inits = list(lambda = lambda,
                                        prob = prob),
                           constants = list(Nmin = Nmin, Nmax = Nmax,
                                            len = len))
          m$calculate()
          MlProbX <- m$getLogProb("x")
          expect_equal(MlProbX, lProbX)

      # Compiled model
          cm <- compileNimble(m)
          cm$calculate()
          CMlProbX <- cm$getLogProb("x")
          expect_equal(CMlProbX, lProbX)

      # Missing values
        xna <- c(1, 0, NA, 3, 2)
        probXna <- dNmixtureAD_s(xna, lambda, prob, Nmin, Nmax, len)
        # Manually calculate the correct answer
        correctProbXna <- 0
        for (N in Nmin:Nmax) {
          correctProbXna <- correctProbXna + dpois(N, lambda) * prod(dbinom(xna, N, prob), na.rm=TRUE)
        }
        expect_equal(probXna, correctProbXna)

        CprobXna <- CdNmixtureAD_s(xna, lambda, prob, Nmin, Nmax, len)
        expect_equal(CprobXna, probXna)

        xna <- c(1,NA,NA,NA,NA)
        correctProbXna <- 0
        for (N in Nmin:Nmax) {
          correctProbXna <- correctProbXna + dpois(N, lambda) * prod(dbinom(xna, N, prob), na.rm=TRUE)
        }
        probXna <- CdNmixtureAD_s(xna, lambda, prob, Nmin, Nmax, len)
        expect_equal(probXna, correctProbXna)

        xna <- as.numeric(rep(NA, 5))
        expect_equal(CdNmixtureAD_s(xna, lambda, prob, Nmin, Nmax, len), 1)
        expect_equal(CdNmixtureAD_s(xna, lambda, prob, Nmin, Nmax, len, log=TRUE), 0)

      # Test imputing value for all NAs
          xNA <- c(NA, NA, NA, NA, NA)
          mNA <- nimbleModel(nc, data = list(x = xNA),
                 inits = list(lambda = lambda,
                              prob = prob),
                 constants = list(Nmin = Nmin, Nmax = Nmax,
                                            len = len))


          mNAConf <- configureMCMC(mNA)
          mNAConf$addMonitors('x')
          mNA_MCMC <- buildMCMC(mNAConf)
          cmNA <- compileNimble(mNA, mNA_MCMC)

          set.seed(0)
          cmNA$mNA_MCMC$run(10)

      # Did the imputed values come back?
          expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x[2]"])))

      # Test simulation code
          nSim <- 10
          xSim <- array(NA, dim = c(nSim, len))
          set.seed(1)
          for (i in 1:nSim) {
            xSim[i,] <- rNmixtureAD_s(1, lambda, prob, Nmin, Nmax, len)
          }

          CrNmixtureAD_s <- compileNimble(rNmixtureAD_s)
          CxSim <- array(NA, dim = c(nSim, len))
          set.seed(1)
          for (i in 1:nSim) {
            CxSim[i,] <- CrNmixtureAD_s(1, lambda, prob, Nmin, Nmax, len)
          }
          expect_identical(xSim, CxSim)

          simNodes <- m$getDependencies(c('prob', 'lambda'), self = FALSE)
          mxSim <- array(NA, dim = c(nSim, len))
          set.seed(1)
          for(i in 1:nSim) {
            m$simulate(simNodes, includeData = TRUE)
            mxSim[i,] <- m$x
          }
          expect_identical(mxSim, xSim)

          CmxSim <- array(NA, dim = c(nSim, len))
          set.seed(1)
          for(i in 1:nSim) {
            cm$simulate(simNodes, includeData = TRUE)
            CmxSim[i,] <- cm$x
          }
          expect_identical(CmxSim, mxSim)
        })



# -----------------------------------------------------------------------------
#### 3. Test dNmixtureAD_BNB_v ####
test_that("dNmixtureAD_BNB_v works",
          {
            # Uncompiled calculation
            x <- c(1, 0, 3, 3, 0)
            lambda <- 8
            theta <- 2
            prob <- c(0.5, 0.3, 0.5, 0.4, 0.1)
            Nmin <- 0
            Nmax <- 250
            len <- 5

            probX <- dNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len)

            # Manually calculate the correct answer
            r <- 1 / theta
            pNB <- 1 / (1 + theta * lambda)
            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dbinom(x, N, prob))
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Dynamic Nmin/Nmax
            expect_error(dynProbX <- dNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin = -1, Nmax = -1, len))

            # Other Nmin / Nmax
            Nmin <- 3
            for(Nmax in 3:6)  {
              dynProbX <- dNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin = Nmin, Nmax = Nmax, len)
              dynCorrectProbX <- 0
              for (N in Nmin:Nmax) {
                dynCorrectProbX <- dynCorrectProbX + dnbinom(N, size = r, prob = pNB) *
                  prod(dbinom(x, N, prob))
              }
              expect_equal(dynProbX, dynCorrectProbX)
            }
            Nmin <- 0
            Nmax <- 250
            # Compilation and compiled calculations
            call_dNmixtureAD_BNB_v <- nimbleFunction(
              name = "t3",
              run=function(x=double(1),
                           lambda=double(),
                           theta=double(),
                           prob=double(1),
                           Nmin = double(0),
                           Nmax = double(0),
                           len=double(),
                           log = integer(0, default=0)) {
                return(dNmixtureAD_BNB_v(x,lambda,theta,prob,Nmin,Nmax,len,log))
                returnType(double())
              }
            )

            # Compilation and compiled calculations
            CdNmixtureAD_BNB_v <- compileNimble(call_dNmixtureAD_BNB_v)
            CprobX <- CdNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len)
            probX <- dNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            lprobX <- dNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            expect_error(CdynProbX <- CdNmixture_BNB_v(x, lambda, theta, prob, Nmin = -1,
                                          Nmax = -1, len))

            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dNmixtureAD_BNB_v(lambda = lambda, prob = prob[1:5],
                                       theta = theta,
                                       Nmin = Nmin, Nmax = Nmax, len = len)

            })

            m <- nimbleModel(code = nc,
                             data = list(x = x),
                             inits = list(lambda = lambda,
                                          prob = prob,
                                          theta = theta),
                             constants = list(Nmin = Nmin, Nmax = Nmax,
                                              len = len))
            m$calculate()
            MlProbX <- m$getLogProb("x")
            expect_equal(MlProbX, lProbX)

            # Compiled model
            cm <- compileNimble(m)
            cm$calculate()
            CMlProbX <- cm$getLogProb("x")
            expect_equal(CMlProbX, lProbX)

            # Missing values
            x <- c(1, 0, NA, 3, 0)
            probX <- dNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len)

            # Manually calculate the correct answer
            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dbinom(x, N, prob), na.rm=TRUE)
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Dynamic Nmin/Nmax
            expect_error(dynProbX <- dNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin = -1, Nmax = -1, len))

            # Other Nmin / Nmax
            Nmin <- 3
            for(Nmax in 3:6)  {
              dynProbX <- dNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin = Nmin, Nmax = Nmax, len)
              dynCorrectProbX <- 0
              for (N in Nmin:Nmax) {
                dynCorrectProbX <- dynCorrectProbX + dnbinom(N, size = r, prob = pNB) *
                  prod(dbinom(x, N, prob), na.rm=TRUE)
              }
              expect_equal(dynProbX, dynCorrectProbX)
            }
            Nmin <- 0
            Nmax <- 250
            # Compilation and compiled calculations
            call_dNmixtureAD_BNB_v <- nimbleFunction(
              name = "t3",
              run=function(x=double(1),
                           lambda=double(),
                           theta=double(),
                           prob=double(1),
                           Nmin = double(0),
                           Nmax = double(0),
                           len=double(),
                           log = integer(0, default=0)) {
                return(dNmixtureAD_BNB_v(x,lambda,theta,prob,Nmin,Nmax,len,log))
                returnType(double())
              }
            )

            # Compilation and compiled calculations
            CdNmixtureAD_BNB_v <- compileNimble(call_dNmixtureAD_BNB_v)
            CprobX <- CdNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len)
            probX <- dNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            lprobX <- dNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            expect_error(CdynProbX <- CdNmixture_BNB_v(x, lambda, theta, prob, Nmin = -1,
                                          Nmax = -1, len))

            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dNmixtureAD_BNB_v(lambda = lambda, prob = prob[1:5],
                                       theta = theta,
                                       Nmin = Nmin, Nmax = Nmax, len = len)

            })

            m <- nimbleModel(code = nc,
                             data = list(x = x),
                             inits = list(lambda = lambda,
                                          prob = prob,
                                          theta = theta),
                             constants = list(Nmin = Nmin, Nmax = Nmax,
                                              len = len))
            m$calculate()
            MlProbX <- m$getLogProb("x")
            expect_equal(MlProbX, lProbX)

            # Compiled model
            cm <- compileNimble(m)
            cm$calculate()
            CMlProbX <- cm$getLogProb("x")
            expect_equal(CMlProbX, lProbX)

            x <- c(1, NA, NA, NA, NA)
            probX <- dNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len)

            # Manually calculate the correct answer
            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dbinom(x, N, prob), na.rm=TRUE)
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            CprobX <- CdNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len)
            probX <- dNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            x <- as.numeric(c(NA, NA, NA, NA, NA))
            probX <- dNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len)

            # Manually calculate the correct answer
            correctProbX <- 1
            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            CprobX <- CdNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len)
            probX <- dNmixtureAD_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            # Test imputing value for all NAs
            xNA <- c(NA, NA, NA, NA, NA)
            mNA <- nimbleModel(nc, data = list(x = xNA),
                               inits = list(lambda = lambda,
                                            theta = theta,
                                            prob = prob),
                               constants = list(Nmin = Nmin, Nmax = Nmax,
                                                len = len))


            mNAConf <- configureMCMC(mNA)
            mNAConf$addMonitors('x')
            mNA_MCMC <- buildMCMC(mNAConf)
            cmNA <- compileNimble(mNA, mNA_MCMC)

            set.seed(0)
            cmNA$mNA_MCMC$run(10)

            # Did the imputed values come back?
            expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x[2]"])))

            # Test simulation code
            nSim <- 10
            xSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for (i in 1:nSim) {
              xSim[i,] <- rNmixtureAD_BNB_v(1, lambda, theta, prob, Nmin, Nmax, len)
            }

            CrNmixtureAD_BNB_v <- compileNimble(rNmixtureAD_BNB_v)
            CxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i,] <- CrNmixtureAD_BNB_v(1, lambda, theta, prob, Nmin, Nmax, len)
            }
            expect_identical(xSim, CxSim)

            simNodes <- m$getDependencies(c('prob', 'lambda'), self = FALSE)
            mxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for(i in 1:nSim) {
              m$simulate(simNodes, includeData = TRUE)
              mxSim[i,] <- m$x
            }
            expect_identical(mxSim, xSim)

            CmxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for(i in 1:nSim) {
              cm$simulate(simNodes, includeData = TRUE)
              CmxSim[i,] <- cm$x
            }
            expect_identical(CmxSim, mxSim)
          })

# -----------------------------------------------------------------------------
#### 4. Test dNmixtureAD_BNB_s ####
test_that("dNmixtureAD_BNB_s works",
          {
            # Uncompiled calculation
            x <- c(1, 0, 3, 3, 0)
            lambda <- 8
            theta <- 2
            prob <- 0.4
            Nmin <- 0
            Nmax <- 250
            len <- 5

            probX <- dNmixtureAD_BNB_s(x, lambda, theta = theta, prob, Nmin, Nmax, len)

            # Manually calculate the correct answer
            r <- 1 / theta
            pNB <- 1 / (1 + theta * lambda)

            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dbinom(x, N, prob))
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Dynamic Nmin/Nmax
            expect_error(dynProbX <- dNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin = -1, Nmax = -1, len))

            # Other Nmin / Nmax
            Nmin <- 3
            for(Nmax in 3:6)  {
              dynProbX <- dNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin = Nmin, Nmax = Nmax, len)
              dynCorrectProbX <- 0
              for (N in Nmin:Nmax) {
                dynCorrectProbX <- dynCorrectProbX + dnbinom(N, size = r, prob = pNB) *
                  prod(dbinom(x, N, prob))
              }
              expect_equal(dynProbX, dynCorrectProbX)
            }
            Nmin <- 0
            Nmax <- 250
            # Compilation and compiled calculations
            call_dNmixtureAD_BNB_s <- nimbleFunction(
              name = "t4",
              run=function(x=double(1),
                           lambda=double(),
                           theta=double(),
                           prob=double(),
                           Nmin = double(0),
                           Nmax = double(0),
                           len=double(),
                           log = integer(0, default=0)) {
                return(dNmixtureAD_BNB_s(x,lambda,theta,prob,Nmin,Nmax,len,log))
                returnType(double())
              }
            )

            # Compilation and compiled calculations
            CdNmixtureAD_BNB_s <- compileNimble(call_dNmixtureAD_BNB_s)
            CprobX <- CdNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            expect_error(CdynProbX <- CdNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin = -1,
                                          Nmax = -1, len))

            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dNmixtureAD_BNB_s(lambda = lambda, prob = prob,
                                       theta = theta,
                                       Nmin = Nmin, Nmax = Nmax, len = len)

            })

            m <- nimbleModel(code = nc,
                             data = list(x = x),
                             inits = list(lambda = lambda,
                                          prob = prob,
                                          theta = theta),
                             constants = list(Nmin = Nmin, Nmax = Nmax,
                                              len = len))
            m$calculate()
            MlProbX <- m$getLogProb("x")
            expect_equal(MlProbX, lProbX)

            # Compiled model
            cm <- compileNimble(m)
            cm$calculate()
            CMlProbX <- cm$getLogProb("x")
            expect_equal(CMlProbX, lProbX)

            # Missing values
            x <- c(1, 0, NA, 3, 0)
            probX <- dNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len)

            # Manually calculate the correct answer
            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dbinom(x, N, prob), na.rm=TRUE)
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Dynamic Nmin/Nmax
            expect_error(dynProbX <- dNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin = -1, Nmax = -1, len))

            # Other Nmin / Nmax
            Nmin <- 3
            for(Nmax in 3:6)  {
              dynProbX <- dNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin = Nmin, Nmax = Nmax, len)
              dynCorrectProbX <- 0
              for (N in Nmin:Nmax) {
                dynCorrectProbX <- dynCorrectProbX + dnbinom(N, size = r, prob = pNB) *
                  prod(dbinom(x, N, prob), na.rm=TRUE)
              }
              expect_equal(dynProbX, dynCorrectProbX)
            }
            Nmin <- 0
            Nmax <- 250
            # Compilation and compiled calculations
            call_dNmixtureAD_BNB_s <- nimbleFunction(
              name = "t3",
              run=function(x=double(1),
                           lambda=double(),
                           theta=double(),
                           prob=double(),
                           Nmin = double(0),
                           Nmax = double(0),
                           len=double(),
                           log = integer(0, default=0)) {
                return(dNmixtureAD_BNB_s(x,lambda,theta,prob,Nmin,Nmax,len,log))
                returnType(double())
              }
            )

            # Compilation and compiled calculations
            CdNmixtureAD_BNB_s <- compileNimble(call_dNmixtureAD_BNB_s)
            CprobX <- CdNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len)
            probX <- dNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            lprobX <- dNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            expect_error(CdynProbX <- CdNmixture_BNB_s(x, lambda, theta, prob, Nmin = -1,
                                          Nmax = -1, len))

            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dNmixtureAD_BNB_s(lambda = lambda, prob = prob,
                                       theta = theta,
                                       Nmin = Nmin, Nmax = Nmax, len = len)

            })

            m <- nimbleModel(code = nc,
                             data = list(x = x),
                             inits = list(lambda = lambda,
                                          prob = prob,
                                          theta = theta),
                             constants = list(Nmin = Nmin, Nmax = Nmax,
                                              len = len))
            m$calculate()
            MlProbX <- m$getLogProb("x")
            expect_equal(MlProbX, lProbX)

            # Compiled model
            cm <- compileNimble(m)
            cm$calculate()
            CMlProbX <- cm$getLogProb("x")
            expect_equal(CMlProbX, lProbX)

            x <- c(1, NA, NA, NA, NA)
            probX <- dNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len)

            # Manually calculate the correct answer
            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dbinom(x, N, prob), na.rm=TRUE)
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            CprobX <- CdNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len)
            probX <- dNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            x <- as.numeric(c(NA, NA, NA, NA, NA))
            probX <- dNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len)

            # Manually calculate the correct answer
            correctProbX <- 1
            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            CprobX <- CdNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len)
            probX <- dNmixtureAD_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            # Test imputing value for all NAs
            xNA <- c(NA, NA, NA, NA, NA)
            mNA <- nimbleModel(nc, data = list(x = xNA),
                               inits = list(lambda = lambda,
                                            theta = theta,
                                            prob = prob),
                               constants = list(Nmin = Nmin, Nmax = Nmax,
                                                len = len))

            mNAConf <- configureMCMC(mNA)
            mNAConf$addMonitors('x')
            mNA_MCMC <- buildMCMC(mNAConf)
            cmNA <- compileNimble(mNA, mNA_MCMC)

            set.seed(0)
            cmNA$mNA_MCMC$run(10)

            # Did the imputed values come back?
            expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x[2]"])))

            # Test simulation code
            nSim <- 10
            xSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for (i in 1:nSim) {
              xSim[i,] <- rNmixtureAD_BNB_s(1, lambda, theta, prob, Nmin, Nmax, len)
            }

            CrNmixtureAD_BNB_s <- compileNimble(rNmixtureAD_BNB_s)
            CxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i,] <- CrNmixtureAD_BNB_s(1, lambda, theta, prob, Nmin, Nmax, len)
            }
            expect_identical(xSim, CxSim)

            simNodes <- m$getDependencies(c('prob', 'lambda'), self = FALSE)
            mxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for(i in 1:nSim) {
              m$simulate(simNodes, includeData = TRUE)
              mxSim[i,] <- m$x
            }
            expect_identical(mxSim, xSim)

            CmxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for(i in 1:nSim) {
              cm$simulate(simNodes, includeData = TRUE)
              CmxSim[i,] <- cm$x
            }
            expect_identical(CmxSim, mxSim)
          })

# -----------------------------------------------------------------------------
#### 5. Test dNmixtureAD_BNB_oneObs ####
test_that("dNmixtureAD_BNB_oneObs works",
          {
            # Uncompiled calculation
            x <- c(1)
            lambda <- 8
            theta <- 2
            prob <- c(0.5)
            Nmin <- 0
            Nmax <- 250
            len <- 5

            probX <- dNmixtureAD_BNB_oneObs(x, lambda, theta, prob, Nmin, Nmax)

            # Manually calculate the correct answer
            r <- 1 / theta
            pNB <- 1 / (1 + theta * lambda)
            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dbinom(x, N, prob))
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BNB_oneObs(x, lambda, theta, prob, Nmin, Nmax, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Dynamic Nmin/Nmax
            expect_error(dynProbX <- dNmixtureAD_BNB_oneObs(x, lambda, theta, prob, Nmin = -1, Nmax = -1))

            # Other Nmin / Nmax
            Nmin <- 3
            for(Nmax in 3:6)  {
              dynProbX <- dNmixtureAD_BNB_oneObs(x, lambda, theta, prob, Nmin = Nmin, Nmax = Nmax)
              dynCorrectProbX <- 0
              for (N in Nmin:Nmax) {
                dynCorrectProbX <- dynCorrectProbX + dnbinom(N, size = r, prob = pNB) *
                  prod(dbinom(x, N, prob))
              }
              expect_equal(dynProbX, dynCorrectProbX)
            }
            Nmin <- 0
            Nmax <- 250

            call_dNmixtureAD_BNB_oneObs <- nimbleFunction(
              name = "t5",
              run=function(x=double(),
                           lambda=double(),
                           theta=double(),
                           prob=double(),
                           Nmin = double(0),
                           Nmax = double(0),
                           log = integer(0, default=0)) {
                return(dNmixtureAD_BNB_oneObs(x,lambda,theta,prob,Nmin,Nmax,log))
                returnType(double())
              }
            )

            # Compilation and compiled calculations
            CdNmixtureAD_BNB_oneObs <- compileNimble(call_dNmixtureAD_BNB_oneObs)
            CprobX <- CdNmixtureAD_BNB_oneObs(x, lambda, theta, prob, Nmin, Nmax)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixtureAD_BNB_oneObs(x, lambda, theta, prob, Nmin, Nmax, log = TRUE)
            expect_equal(ClProbX, lProbX)

            expect_error(CdynProbX <- CdNmixtureAD_BNB_oneObs(x, lambda, theta, prob, Nmin = -1, Nmax = -1))

            # Use in Nimble model
            nc <- nimbleCode({
              x ~ dNmixtureAD_BNB_oneObs(lambda = lambda, prob = prob,
                                       theta = theta,
                                       Nmin = Nmin, Nmax = Nmax)
            })

            m <- nimbleModel(code = nc,
                             data = list(x = x),
                             inits = list(lambda = lambda,
                                          prob = prob,
                                          theta = theta),
                             constants = list(Nmin = Nmin, Nmax = Nmax))
            m$calculate()
            MlProbX <- m$getLogProb("x")
            expect_equal(MlProbX, lProbX)

            # Compiled model
            cm <- compileNimble(m)
            cm$calculate()
            CMlProbX <- cm$getLogProb("x")
            expect_equal(CMlProbX, lProbX)

            # Missing values
            x <- as.numeric(NA)
            probX <- dNmixtureAD_BNB_oneObs(x, lambda, theta, prob, Nmin, Nmax)
            expect_equal(probX, 1)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BNB_oneObs(x, lambda, theta, prob, Nmin, Nmax, log = TRUE)
            lCorrectProbX <- log(1)
            expect_equal(lProbX, lCorrectProbX)

            # Compilation and compiled calculations
            CprobX <- CdNmixtureAD_BNB_oneObs(x, lambda, theta, prob, Nmin, Nmax)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixtureAD_BNB_oneObs(x, lambda, theta, prob, Nmin, Nmax, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Test imputing value for all NAs
            xNA <- NA
            mNA <- nimbleModel(nc, data = list(x = xNA),
                               inits = list(lambda = lambda,
                                            theta = theta,
                                            prob = prob),
                               constants = list(Nmin = Nmin, Nmax = Nmax))

            mNAConf <- configureMCMC(mNA)
            mNAConf$addMonitors('x')
            mNA_MCMC <- buildMCMC(mNAConf)
            cmNA <- compileNimble(mNA, mNA_MCMC)

            set.seed(0)
            cmNA$mNA_MCMC$run(10)

            # Did the imputed values come back?
            expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x"])))

            # Test simulation code
            nSim <- 10
            xSim <- numeric(nSim)
            set.seed(1)
            for (i in 1:nSim) {
              xSim[i] <- rNmixtureAD_BNB_oneObs(1, lambda, theta, prob, Nmin, Nmax)
            }

            CrNmixtureAD_BNB_oneObs <- compileNimble(rNmixtureAD_BNB_oneObs)
            CxSim <- numeric(nSim)
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i] <- CrNmixtureAD_BNB_oneObs(1, lambda, theta, prob, Nmin, Nmax)
            }
            expect_identical(xSim, CxSim)

            simNodes <- m$getDependencies(c('prob', 'lambda'), self = FALSE)
            mxSim <- numeric(nSim)
            set.seed(1)
            for(i in 1:nSim) {
              m$simulate(simNodes, includeData = TRUE)
              mxSim[i] <- m$x
            }
            expect_identical(mxSim, xSim)

            CmxSim <- numeric(nSim)
            set.seed(1)
            for(i in 1:nSim) {
              cm$simulate(simNodes, includeData = TRUE)
              CmxSim[i] <- cm$x
            }
            expect_identical(CmxSim, mxSim)
          })


# -----------------------------------------------------------------------------
#### 6. Test dNmixtureAD_BBP_v ####
test_that("dNmixtureAD_BBP_v works",
          {
            # Uncompiled calculation
            x <- c(1, 0, 3, 3, 0)
            lambda <- 8
            s <- 2
            prob <- c(0.5, 0.3, 0.5, 0.4, 0.1)
            Nmin <- max(x)
            Nmax <- 250
            len <- 5

            probX <- dNmixtureAD_BBP_v(x, lambda, prob, s, Nmin, Nmax, len)

            # Manually calculate the correct answer
            alpha <- prob * s
            beta <- s - prob * s

            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dpois(N, lambda) *
                prod(dBetaBinom_v(x, N, shape1 = alpha, shape2 = beta))
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BBP_v(x, lambda, prob, s, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Other Nmin / Nmax
            Nmin <- 3
            for(Nmax in 3:6)  {
              dynProbX <- dNmixtureAD_BBP_v(x, lambda, prob, s, Nmin = Nmin, Nmax = Nmax, len)
              dynCorrectProbX <- 0
              for (N in Nmin:Nmax) {
                dynCorrectProbX <- dynCorrectProbX + dpois(N, lambda) *
                  prod(dBetaBinom_v(x, N, shape1 = alpha, shape2 = beta))
              }
              expect_equal(dynProbX, dynCorrectProbX)
            }
            Nmin <- 0
            Nmax <- 250
            # Compilation and compiled calculations
            call_dNmixtureAD_BBP_v <- nimbleFunction(
              name = "t6",
              run=function(x=double(1),
                           lambda=double(),
                           prob=double(1),
                           s=double(),
                           Nmin = double(0),
                           Nmax = double(0),
                           len=double(),
                           log = integer(0, default=0)) {
                return(dNmixtureAD_BBP_v(x,lambda,prob,s,Nmin,Nmax,len,log))
                returnType(double())
              }
            )
            # Compilation and compiled calculations
            CdNmixtureAD_BBP_v <- compileNimble(call_dNmixtureAD_BBP_v)
            CprobX <- CdNmixtureAD_BBP_v(x, lambda, prob, s, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixtureAD_BBP_v(x, lambda, prob, s, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Dynamic Nmin / Nmax isn't allowed
            expect_error({
              dNmixtureAD_BBP_v(x, lambda, prob, s, Nmin = -1, Nmax = -1, len)
            })
            expect_error({
              dNmixtureAD_BBP_v(x, lambda, prob, s, Nmin = -1, Nmax, len)
            })
            expect_error({
              dNmixtureAD_BBP_v(x, lambda, prob, s, Nmin, Nmax = -1, len)
            })
            expect_error({
              CdNmixtureAD_BBP_v(x, lambda, prob, s, Nmin = -1, Nmax = -1, len)
            })
            expect_error({
              CdNmixtureAD_BBP_v(x, lambda, prob, s, Nmin = -1, Nmax, len)
            })
            expect_error({
              CdNmixtureAD_BBP_v(x, lambda, prob, s, Nmin, Nmax = -1, len)
            })


            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dNmixtureAD_BBP_v(lambda = lambda, prob = prob[1:5],
                                       s = s,
                                       Nmin = Nmin, Nmax = Nmax, len = len)
            })

            m <- nimbleModel(code = nc,
                             data = list(x = x),
                             inits = list(lambda = lambda,
                                          prob = prob,
                                          s = s),
                             constants = list(Nmin = Nmin, Nmax = Nmax,
                                              len = len))
            m$calculate()
            MlProbX <- m$getLogProb("x")
            expect_equal(MlProbX, lProbX)

            # Compiled model
            cm <- compileNimble(m)
            cm$calculate()
            CMlProbX <- cm$getLogProb("x")
            expect_equal(CMlProbX, lProbX)

            # Missing values
            # Uncompiled calculation
            x <- c(1, 0, NA, 3, 0)
            Nmin <- max(x, na.rm=TRUE)
            Nmax <- 250
            probX <- dNmixtureAD_BBP_v(x, lambda, prob, s, Nmin, Nmax, len)

            # Manually calculate the correct answer
            alpha <- prob * s
            beta <- s - prob * s

            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dpois(N, lambda) *
                prod(dBetaBinom_v(x, N, shape1 = alpha, shape2 = beta), na.rm=TRUE)
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BBP_v(x, lambda, prob, s, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Other Nmin / Nmax
            Nmin <- 3
            for(Nmax in 3:6)  {
              dynProbX <- dNmixtureAD_BBP_v(x, lambda, prob, s, Nmin = Nmin, Nmax = Nmax, len)
              dynCorrectProbX <- 0
              for (N in Nmin:Nmax) {
                dynCorrectProbX <- dynCorrectProbX + dpois(N, lambda) *
                  prod(dBetaBinom_v(x, N, shape1 = alpha, shape2 = beta))
              }
              expect_equal(dynProbX, dynCorrectProbX)
            }
            Nmin <- 0
            Nmax <- 250
            # Compilation and compiled calculations
            CprobX <- CdNmixtureAD_BBP_v(x, lambda, prob, s, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixtureAD_BBP_v(x, lambda, prob, s, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            m <- nimbleModel(code = nc,
                             data = list(x = x),
                             inits = list(lambda = lambda,
                                          prob = prob,
                                          s = s),
                             constants = list(Nmin = Nmin, Nmax = Nmax,
                                              len = len))
            m$calculate()
            MlProbX <- m$getLogProb("x")
            expect_equal(MlProbX, lProbX)

            # Compiled model
            cm <- compileNimble(m)
            cm$calculate()
            CMlProbX <- cm$getLogProb("x")
            expect_equal(CMlProbX, lProbX)

            x <- as.numeric(rep(NA,5))
            probX <- dNmixtureAD_BBP_v(x, lambda, prob, s, Nmin, Nmax, len)
            expect_equal(probX, 1)

            # Compilation and compiled calculations
            CprobX <- CdNmixtureAD_BBP_v(x, lambda, prob, s, Nmin, Nmax, len)
            expect_equal(CprobX, 1)

            # Test imputing value for all NAs
            xNA <- c(NA, NA, NA, NA, NA)
            mNA <- nimbleModel(nc, data = list(x = xNA),
                               inits = list(lambda = lambda,
                                            s = s,
                                            prob = prob),
                               constants = list(Nmin = Nmin, Nmax = Nmax,
                                                len = len))

            mNAConf <- configureMCMC(mNA)
            mNAConf$addMonitors('x')
            mNA_MCMC <- buildMCMC(mNAConf)
            cmNA <- compileNimble(mNA, mNA_MCMC)

            set.seed(0)
            cmNA$mNA_MCMC$run(10)

            # Did the imputed values come back?
            expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x[2]"])))

            # Test simulation code
            nSim <- 10
            xSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for (i in 1:nSim) {
              xSim[i,] <- rNmixtureAD_BBP_v(1, lambda, prob, s, Nmin, Nmax, len)
            }

            CrNmixtureAD_BBP_v <- compileNimble(rNmixtureAD_BBP_v)
            CxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i,] <- CrNmixtureAD_BBP_v(1, lambda, prob, s, Nmin, Nmax, len)
            }
            expect_equal(xSim, CxSim)

            simNodes <- m$getDependencies(c('prob', 'lambda'), self = FALSE)
            mxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for(i in 1:nSim) {
              m$simulate(simNodes, includeData = TRUE)
              mxSim[i,] <- m$x
            }
            expect_equal(mxSim, xSim)

            CmxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for(i in 1:nSim) {
              cm$simulate(simNodes, includeData = TRUE)
              CmxSim[i,] <- cm$x
            }
            expect_identical(CmxSim, mxSim)
          })

# -----------------------------------------------------------------------------
#### 7. Test dNmixtureAD_BBP_s ####
test_that("dNmixtureAD_BBP_s works",
          {
            # Uncompiled calculation
            x <- c(1, 0, 3, 3, 0)
            lambda <- 8
            s <- 2
            prob <- 0.4
            Nmin <- max(x)
            Nmax <- 250
            len <- 5

            probX <- dNmixtureAD_BBP_s(x, lambda, prob, s, Nmin, Nmax, len)

            # Manually calculate the correct answer
            alpha <- prob * s
            beta <- s - prob * s

            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dpois(N, lambda) *
                prod(dBetaBinom_s(x, N,
                                shape1 = alpha, shape2 = beta))
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BBP_s(x, lambda, prob, s, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Other Nmin / Nmax
            Nmin <- 3
            for(Nmax in 3:6)  {
              dynProbX <- dNmixtureAD_BBP_s(x, lambda, prob, s, Nmin = Nmin, Nmax = Nmax, len)
              dynCorrectProbX <- 0
              for (N in Nmin:Nmax) {
                dynCorrectProbX <- dynCorrectProbX + dpois(N, lambda) *
                  prod(dBetaBinom_s(x, N, shape1 = alpha, shape2 = beta))
              }
              expect_equal(dynProbX, dynCorrectProbX)
            }
            Nmin <- 0
            Nmax <- 250

            # Compilation and compiled calculations
            call_dNmixtureAD_BBP_s <- nimbleFunction(
              name = "t7",
              run=function(x=double(1),
                           lambda=double(),
                           prob=double(),
                           s=double(),
                           Nmin = double(0),
                           Nmax = double(0),
                           len=double(),
                           log = integer(0, default=0)) {
                return(dNmixtureAD_BBP_s(x,lambda,prob,s,Nmin,Nmax,len,log))
                returnType(double())
              }
            )
            # Compilation and compiled calculations
            CdNmixtureAD_BBP_s <- compileNimble(call_dNmixtureAD_BBP_s)
            CprobX <- CdNmixtureAD_BBP_s(x, lambda, prob, s, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixtureAD_BBP_s(x, lambda, prob, s, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)


            # Dynamic Nmin / Nmax isn't allowed
            expect_error({
              dNmixtureAD_BBP_s(x, lambda, prob, s, Nmin = -1, Nmax = -1, len)
            })
            expect_error({
              dNmixtureAD_BBP_s(x, lambda, prob, s, Nmin = -1, Nmax, len)
            })
            expect_error({
              dNmixtureAD_BBP_s(x, lambda, prob, s, Nmin, Nmax = -1, len)
            })
            expect_error({
              CdNmixtureAD_BBP_s(x, lambda, prob, s, Nmin = -1, Nmax = -1, len)
            })
            expect_error({
              CdNmixtureAD_BBP_s(x, lambda, prob, s, Nmin = -1, Nmax, len)
            })
            expect_error({
              CdNmixtureAD_BBP_s(x, lambda, prob, s, Nmin, Nmax = -1, len)
            })


            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dNmixtureAD_BBP_s(lambda = lambda, prob = prob,
                                       s = s,
                                       Nmin = Nmin, Nmax = Nmax, len = len)

            })

            m <- nimbleModel(code = nc,
                             data = list(x = x),
                             inits = list(lambda = lambda,
                                          prob = prob,
                                          s = s),
                             constants = list(Nmin = Nmin, Nmax = Nmax,
                                              len = len))
            m$calculate()
            MlProbX <- m$getLogProb("x")
            expect_equal(MlProbX, lProbX)

            # Compiled model
            cm <- compileNimble(m)
            cm$calculate()
            CMlProbX <- cm$getLogProb("x")
            expect_equal(CMlProbX, lProbX)

            # Missing values
            # Uncompiled calculation
            x <- c(1, 0, NA, 3, 0)
            Nmin <- max(x, na.rm=TRUE)
            
            probX <- dNmixtureAD_BBP_s(x, lambda, prob, s, Nmin, Nmax, len)

            # Manually calculate the correct answer
            alpha <- prob * s
            beta <- s - prob * s

            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dpois(N, lambda) *
                prod(dBetaBinom_s(x, N,
                                shape1 = alpha, shape2 = beta), na.rm=TRUE)
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BBP_s(x, lambda, prob, s, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Other Nmin / Nmax
            Nmin <- 3
            for(Nmax in 3:6)  {
              dynProbX <- dNmixtureAD_BBP_s(x, lambda, prob, s, Nmin = Nmin, Nmax = Nmax, len)
              dynCorrectProbX <- 0
              for (N in Nmin:Nmax) {
                dynCorrectProbX <- dynCorrectProbX + dpois(N, lambda) *
                  prod(dBetaBinom_s(x, N, shape1 = alpha, shape2 = beta))
              }
              expect_equal(dynProbX, dynCorrectProbX)
            }
            Nmin <- 0
            Nmax <- 250

            # Compilation and compiled calculations
            CprobX <- CdNmixtureAD_BBP_s(x, lambda, prob, s, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixtureAD_BBP_s(x, lambda, prob, s, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            m <- nimbleModel(code = nc,
                             data = list(x = x),
                             inits = list(lambda = lambda,
                                          prob = prob,
                                          s = s),
                             constants = list(Nmin = Nmin, Nmax = Nmax,
                                              len = len))
            m$calculate()
            MlProbX <- m$getLogProb("x")
            expect_equal(MlProbX, lProbX)

            # Compiled model
            cm <- compileNimble(m)
            cm$calculate()
            CMlProbX <- cm$getLogProb("x")
            expect_equal(CMlProbX, lProbX)

            x <- as.numeric(rep(NA, 5))
            probX <- dNmixtureAD_BBP_s(x, lambda, prob, s, Nmin, Nmax, len)
            expect_equal(probX, 1)

              # Compilation and compiled calculations
            CprobX <- CdNmixtureAD_BBP_s(x, lambda, prob, s, Nmin, Nmax, len)
            expect_equal(CprobX, 1)

            # Test imputing value for all NAs
            xNA <- c(NA, NA, NA, NA, NA)
            mNA <- nimbleModel(nc, data = list(x = xNA),
                               inits = list(lambda = lambda,
                                            s = s,
                                            prob = prob),
                               constants = list(Nmin = Nmin, Nmax = Nmax,
                                                len = len))


            mNAConf <- configureMCMC(mNA)
            mNAConf$addMonitors('x')
            mNA_MCMC <- buildMCMC(mNAConf)
            cmNA <- compileNimble(mNA, mNA_MCMC)

            set.seed(0)
            cmNA$mNA_MCMC$run(10)

            # Did the imputed values come back?
            expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x[2]"])))

            # Test simulation code
            nSim <- 10
            xSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for (i in 1:nSim) {
              xSim[i,] <- rNmixtureAD_BBP_s(1, lambda, prob, s, Nmin, Nmax, len)
            }

            CrNmixtureAD_BBP_s <- compileNimble(rNmixtureAD_BBP_s)
            CxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i,] <- CrNmixtureAD_BBP_s(1, lambda, prob, s, Nmin, Nmax, len)
            }
            expect_equal(xSim, CxSim)

            simNodes <- m$getDependencies(c('prob', 'lambda'), self = FALSE)
            mxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for(i in 1:nSim) {
              m$simulate(simNodes, includeData = TRUE)
              mxSim[i,] <- m$x
            }
            expect_equal(mxSim, xSim)

            CmxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for(i in 1:nSim) {
              cm$simulate(simNodes, includeData = TRUE)
              CmxSim[i,] <- cm$x
            }
            expect_identical(CmxSim, mxSim)
          })


# -----------------------------------------------------------------------------
#### 8. Test dNmixtureAD_BBP_oneObs ####
test_that("dNmixtureAD_BBP_oneObs works",
          {
            # Uncompiled calculation
            x <- c(1)
            lambda <- 8
            s <- 2
            prob <- c(0.5)
            Nmin <- max(x)
            Nmax <- 250
            len <- 5

            probX <- dNmixtureAD_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax)

            # Manually calculate the correct answer
            alpha <- prob * s
            beta <- s - prob * s
            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dpois(N, lambda) *
                prod(dBetaBinom_s(x, N, alpha, beta))
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Other Nmin / Nmax
            Nmin <- 3
            for(Nmax in 3:6)  {
              dynProbX <- dNmixtureAD_BBP_oneObs(x, lambda, prob, s, Nmin = Nmin, Nmax = Nmax)
              dynCorrectProbX <- 0
              for (N in Nmin:Nmax) {
                dynCorrectProbX <- dynCorrectProbX + dpois(N, lambda) *
                  prod(dBetaBinom_s(x, N, shape1 = alpha, shape2 = beta))
              }
              expect_equal(dynProbX, dynCorrectProbX)
            }
            Nmin <- 0
            Nmax <- 250

            # Compilation and compiled calculations
            call_dNmixtureAD_BBP_oneObs <- nimbleFunction(
              name = "t8",
              run=function(x=double(),
                           lambda=double(),
                           prob=double(),
                           s=double(),
                           Nmin = double(0),
                           Nmax = double(0),
                           log = integer(0, default=0)) {
                return(dNmixtureAD_BBP_oneObs(x,lambda,prob,s,Nmin,Nmax,log))
                returnType(double())
              }
            )

            # Compilation and compiled calculations
            CdNmixtureAD_BBP_oneObs <- compileNimble(call_dNmixtureAD_BBP_oneObs)
            CprobX <- CdNmixtureAD_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixtureAD_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Dynamic Nmin / Nmax isn't allowed
            expect_error({
              dNmixtureAD_BBP_oneObs(x, lambda, prob, s, Nmin = -1, Nmax = -1)
            })
            expect_error({
              dNmixtureAD_BBP_oneObs(x, lambda, prob, s, Nmin = -1, Nmax)
            })
            expect_error({
              dNmixtureAD_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax = -1)
            })
            expect_error({
              CdNmixtureAD_BBP_oneObs(x, lambda, prob, s, Nmin = -1, Nmax = -1)
            })
            expect_error({
              CdNmixtureAD_BBP_oneObs(x, lambda, prob, s, Nmin = -1, Nmax)
            })
            expect_error({
              CdNmixtureAD_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax = -1)
            })

            # Use in Nimble model
            nc <- nimbleCode({
              x ~ dNmixtureAD_BBP_oneObs(lambda = lambda, prob = prob,
                                  s = s,
                                  Nmin = Nmin, Nmax = Nmax)

            })

            m <- nimbleModel(code = nc,
                             data = list(x = x),
                             inits = list(lambda = lambda,
                                          prob = prob,
                                          s = s),
                             constants = list(Nmin = Nmin, Nmax = Nmax))
            m$calculate()
            MlProbX <- m$getLogProb("x")
            expect_equal(MlProbX, lProbX)

            # Compiled model
            cm <- compileNimble(m)
            cm$calculate()
            CMlProbX <- cm$getLogProb("x")
            expect_equal(CMlProbX, lProbX)

            # Missing value
            # Uncompiled calculation
            x <- as.numeric(NA)
            lambda <- 8
            s <- 2
            prob <- c(0.5)
            Nmin <- 0
            Nmax <- 250
            len <- 5

            probX <- dNmixtureAD_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax)

            expect_equal(probX, 1)

            # Compilation and compiled calculations
            CprobX <- CdNmixtureAD_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax)
            expect_equal(CprobX, 1)

            # Test imputing value for all NAs
            xNA <- NA
            mNA <- nimbleModel(nc, data = list(x = xNA),
                               inits = list(lambda = lambda,
                                            s = s,
                                            prob = prob),
                               constants = list(Nmin = Nmin, Nmax = Nmax))


            mNAConf <- configureMCMC(mNA)
            mNAConf$addMonitors('x')
            mNA_MCMC <- buildMCMC(mNAConf)
            cmNA <- compileNimble(mNA, mNA_MCMC)

            set.seed(0)
            cmNA$mNA_MCMC$run(10)

            # Did the imputed values come back?
            expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x"])))

            # Test simulation code
            nSim <- 10
            xSim <- numeric(nSim)
            set.seed(1)
            for (i in 1:nSim) {
              xSim[i] <- rNmixtureAD_BBP_oneObs(1, lambda, prob, s, Nmin, Nmax)
            }

            CrNmixtureAD_BBP_oneObs <- compileNimble(rNmixtureAD_BBP_oneObs)
            CxSim <- numeric(nSim)
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i] <- CrNmixtureAD_BBP_oneObs(1, lambda, prob, s, Nmin, Nmax)
            }
            expect_identical(xSim, CxSim)

            simNodes <- m$getDependencies(c('prob', 'lambda'), self = FALSE)
            mxSim <- numeric(nSim)
            set.seed(1)
            for(i in 1:nSim) {
              m$simulate(simNodes, includeData = TRUE)
              mxSim[i] <- m$x
            }
            expect_identical(mxSim, xSim)

            CmxSim <- numeric(nSim)
            set.seed(1)
            for(i in 1:nSim) {
              cm$simulate(simNodes, includeData = TRUE)
              CmxSim[i] <- cm$x
            }
            expect_identical(CmxSim, mxSim)
          })




# -----------------------------------------------------------------------------
#### 9. Test dNmixtureAD_BBNB_v ####
test_that("dNmixtureAD_BBNB_v works",
          {
            # Uncompiled calculation
            x <- c(1, 0, 3, 3, 0)
            lambda <- 8
            s <- 2
            theta <- 1.5
            prob <- c(0.5, 0.3, 0.5, 0.4, 0.1)
            Nmin <- max(x)
            Nmax <- 250
            len <- 5

            probX <- dNmixtureAD_BBNB_v(x, lambda, theta, prob, s, Nmin, Nmax, len)

            # Manually calculate the correct answer
            alpha <- prob * s
            beta <- s - prob * s
            r <- 1 / theta
            pNB <- 1 / (1 + theta * lambda)

            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dBetaBinom_v(x, N, shape1 = alpha, shape2 = beta))
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BBNB_v(x, lambda, theta, prob, s, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Other Nmin / Nmax
            Nmin <- 3
            for(Nmax in 3:6)  {
              dynProbX <- dNmixtureAD_BBNB_v(x, lambda, theta, prob, s, Nmin = Nmin, Nmax = Nmax, len)
              dynCorrectProbX <- 0
              for (N in Nmin:Nmax) {
                dynCorrectProbX <- dynCorrectProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dBetaBinom_v(x, N, shape1 = alpha, shape2 = beta))
              }
              expect_equal(dynProbX, dynCorrectProbX)
            }
            Nmin <- 0
            Nmax <- 250

            # Compilation and compiled calculations
            call_dNmixtureAD_BBNB_v <- nimbleFunction(
              name = "t9",
              run=function(x=double(1),
                           lambda=double(),
                           theta=double(),
                           prob=double(1),
                           s=double(),
                           Nmin = double(0),
                           Nmax = double(0),
                           len=double(),
                           log = integer(0, default=0)) {
                return(dNmixtureAD_BBNB_v(x,lambda,theta,prob,s,Nmin,Nmax,len,log))
                returnType(double())
              }
            )
            # Compilation and compiled calculations
            CdNmixtureAD_BBNB_v <- compileNimble(call_dNmixtureAD_BBNB_v)
            CprobX <- CdNmixtureAD_BBNB_v(x, lambda, theta, prob, s, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixtureAD_BBNB_v(x, lambda, theta, prob, s, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Dynamic Nmin / Nmax isn't allowed
            expect_error({
              dNmixtureAD_BBNB_v(x, lambda, theta, prob, s, Nmin = -1, Nmax = -1, len)
            })
            expect_error({
              dNmixtureAD_BBNB_v(x, lambda, theta, prob, s, Nmin, Nmax = -1, len)
            })
            expect_error({
              dNmixtureAD_BBNB_v(x, lambda, theta, prob, s, Nmin = -1, Nmax, len)
            })
            expect_error({
              CdNmixtureAD_BBNB_v(x, lambda, theta, prob, s, Nmin = -1, Nmax = -1, len)
            })
            expect_error({
              CdNmixtureAD_BBNB_v(x, lambda, theta, prob, s, Nmin, Nmax = -1, len)
            })
            expect_error({
              CdNmixtureAD_BBNB_v(x, lambda, theta, prob, s, Nmin = -1, Nmax, len)
            })

            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dNmixtureAD_BBNB_v(lambda = lambda, prob = prob[1:5],
                                       theta = theta, s = s,
                                       Nmin = Nmin, Nmax = Nmax, len = len)

            })

            m <- nimbleModel(code = nc,
                             data = list(x = x),
                             inits = list(lambda = lambda,
                                          prob = prob,
                                          s = s, theta = theta),
                             constants = list(Nmin = Nmin, Nmax = Nmax,
                                              len = len))
            m$calculate()
            MlProbX <- m$getLogProb("x")
            expect_equal(MlProbX, lProbX)

            # Compiled model
            cm <- compileNimble(m)
            cm$calculate()
            CMlProbX <- cm$getLogProb("x")
            expect_equal(CMlProbX, lProbX)

            # Missing values
            # Uncompiled calculation
            x <- c(1, 0, NA, 3, 0)
            Nmin <- max(x, na.rm=TRUE)
            probX <- dNmixtureAD_BBNB_v(x, lambda, theta, prob, s, Nmin, Nmax, len)

            # Manually calculate the correct answer
            alpha <- prob * s
            beta <- s - prob * s
            r <- 1 / theta
            pNB <- 1 / (1 + theta * lambda)

            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dBetaBinom_v(x, N, shape1 = alpha, shape2 = beta), na.rm=TRUE)
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BBNB_v(x, lambda, theta, prob, s, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Other Nmin / Nmax
            Nmin <- 3
            for(Nmax in 3:6)  {
              dynProbX <- dNmixtureAD_BBNB_v(x, lambda, theta, prob, s, Nmin = Nmin, Nmax = Nmax, len)
              dynCorrectProbX <- 0
              for (N in Nmin:Nmax) {
                dynCorrectProbX <- dynCorrectProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dBetaBinom_v(x, N, shape1 = alpha, shape2 = beta))
              }
              expect_equal(dynProbX, dynCorrectProbX)
            }
            Nmin <- 0
            Nmax <- 250

            CprobX <- CdNmixtureAD_BBNB_v(x, lambda, theta, prob, s, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixtureAD_BBNB_v(x, lambda, theta, prob, s, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dNmixtureAD_BBNB_v(lambda = lambda, prob = prob[1:5],
                                       theta = theta, s = s,
                                       Nmin = Nmin, Nmax = Nmax, len = len)

            })

            m <- nimbleModel(code = nc,
                             data = list(x = x),
                             inits = list(lambda = lambda,
                                          prob = prob,
                                          s = s, theta = theta),
                             constants = list(Nmin = Nmin, Nmax = Nmax,
                                              len = len))
            m$calculate()
            MlProbX <- m$getLogProb("x")
            expect_equal(MlProbX, lProbX)

            # Compiled model
            cm <- compileNimble(m)
            cm$calculate()
            CMlProbX <- cm$getLogProb("x")
            expect_equal(CMlProbX, lProbX)

            # All missing
            x <- as.numeric(rep(NA, 5))
            probX <- dNmixtureAD_BBNB_v(x, lambda, theta, prob, s, Nmin, Nmax, len)
            expect_equal(probX, 1)

            CprobX <- CdNmixtureAD_BBNB_v(x, lambda, theta, prob, s, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            # Test imputing value for all NAs
            xNA <- c(NA, NA, NA, NA, NA)
            mNA <- nimbleModel(nc, data = list(x = xNA),
                               inits = list(lambda = lambda,
                                            s = s, theta = theta,
                                            prob = prob),
                               constants = list(Nmin = Nmin, Nmax = Nmax,
                                                len = len))

            mNAConf <- configureMCMC(mNA)
            mNAConf$addMonitors('x')
            mNA_MCMC <- buildMCMC(mNAConf)
            cmNA <- compileNimble(mNA, mNA_MCMC)

            set.seed(0)
            cmNA$mNA_MCMC$run(10)

            # Did the imputed values come back?
            expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x[2]"])))

            # Test simulation code
            nSim <- 10
            xSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for (i in 1:nSim) {
              xSim[i,] <- rNmixtureAD_BBNB_v(1, lambda, theta, prob, s, Nmin, Nmax, len)
            }

            CrNmixtureAD_BBNB_v <- compileNimble(rNmixtureAD_BBNB_v)
            CxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i,] <- CrNmixtureAD_BBNB_v(1, lambda, theta, prob, s, Nmin, Nmax, len)
            }
            expect_equal(xSim, CxSim)

            simNodes <- m$getDependencies(c('prob', 'lambda'), self = FALSE)
            mxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for(i in 1:nSim) {
              m$simulate(simNodes, includeData = TRUE)
              mxSim[i,] <- m$x
            }
            expect_equal(mxSim, xSim)

            CmxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for(i in 1:nSim) {
              cm$simulate(simNodes, includeData = TRUE)
              CmxSim[i,] <- cm$x
            }
            expect_identical(CmxSim, mxSim)
          })

# -----------------------------------------------------------------------------
#### 10. Test dNmixtureAD_BBNB_s ####
test_that("dNmixtureAD_BBNB_s works",
          {
            # Uncompiled calculation
            x <- c(1, 0, 3, 3, 0)
            lambda <- 8
            s <- 2
            theta <- 1.5
            prob <- 0.4
            Nmin <- max(x)
            Nmax <- 250
            len <- 5

            probX <- dNmixtureAD_BBNB_s(x, lambda, theta, prob, s, Nmin, Nmax, len)

            # Manually calculate the correct answer
            alpha <- prob * s
            beta <- s - prob * s
            r <- 1 / theta
            pNB <- 1 / (1 + theta * lambda)

            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dBetaBinom_s(x, N, shape1 = alpha, shape2 = beta))
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BBNB_s(x, lambda, theta, prob, s, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)
            # Other Nmin / Nmax
            Nmin <- 3
            for(Nmax in 3:6)  {
              dynProbX <- dNmixtureAD_BBNB_s(x, lambda, theta, prob, s, Nmin = Nmin, Nmax = Nmax, len)
              dynCorrectProbX <- 0
              for (N in Nmin:Nmax) {
                dynCorrectProbX <- dynCorrectProbX + dnbinom(N, size = r, prob = pNB) *
                  prod(dBetaBinom_s(x, N, shape1 = alpha, shape2 = beta))
              }
              expect_equal(dynProbX, dynCorrectProbX)
            }
            Nmin <- 0
            Nmax <- 250

            # Compilation and compiled calculations
            call_dNmixtureAD_BBNB_s <- nimbleFunction(
              name = "t10",
              run=function(x=double(1),
                           lambda=double(),
                           theta=double(),
                           prob=double(),
                           s=double(),
                           Nmin = double(0),
                           Nmax = double(0),
                           len=double(),
                           log = integer(0, default=0)) {
                return(dNmixtureAD_BBNB_s(x,lambda,theta,prob,s,Nmin,Nmax,len,log))
                returnType(double())
              }
            )
            # Compilation and compiled calculations
            CdNmixtureAD_BBNB_s <- compileNimble(call_dNmixtureAD_BBNB_s)
            CprobX <- CdNmixtureAD_BBNB_s(x, lambda, theta, prob, s, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixtureAD_BBNB_s(x, lambda, theta, prob, s, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Dynamic Nmin / Nmax isn't allowed
            expect_error({
              dNmixtureAD_BBNB_s(x, lambda, theta, prob, s, Nmin = -1, Nmax = -1, len)
            })
            expect_error({
              dNmixtureAD_BBNB_s(x, lambda, theta, prob, s, Nmin, Nmax = -1, len)
            })
            expect_error({
              dNmixtureAD_BBNB_s(x, lambda, theta, prob, s, Nmin = -1, Nmax, len)
            })
            expect_error({
              CdNmixtureAD_BBNB_s(x, lambda, theta, prob, s, Nmin = -1, Nmax = -1, len)
            })
            expect_error({
              CdNmixtureAD_BBNB_s(x, lambda, theta, prob, s, Nmin, Nmax = -1, len)
            })
            expect_error({
              CdNmixtureAD_BBNB_s(x, lambda, theta, prob, s, Nmin = -1, Nmax, len)
            })
            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dNmixtureAD_BBNB_s(lambda = lambda, prob = prob,
                                        theta = theta, s = s,
                                        Nmin = Nmin, Nmax = Nmax, len = len)
            })

            m <- nimbleModel(code = nc,
                             data = list(x = x),
                             inits = list(lambda = lambda,
                                          prob = prob,
                                          s = s, theta = theta),
                             constants = list(Nmin = Nmin, Nmax = Nmax,
                                              len = len))
            m$calculate()
            MlProbX <- m$getLogProb("x")
            expect_equal(MlProbX, lProbX)

            # Compiled model
            cm <- compileNimble(m)
            cm$calculate()
            CMlProbX <- cm$getLogProb("x")
            expect_equal(CMlProbX, lProbX)

            # Missing values
            # Uncompiled calculation
            x <- c(1, 0, NA, 3, 0)
            Nmin <- max(x, na.rm=TRUE)

            probX <- dNmixtureAD_BBNB_s(x, lambda, theta, prob, s, Nmin, Nmax, len)

            # Manually calculate the correct answer
            alpha <- prob * s
            beta <- s - prob * s
            r <- 1 / theta
            pNB <- 1 / (1 + theta * lambda)

            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dBetaBinom_s(x, N, shape1 = alpha, shape2 = beta), na.rm=TRUE)
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BBNB_s(x, lambda, theta, prob, s, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)
            # Other Nmin / Nmax
            Nmin <- 3
            for(Nmax in 3:6)  {
              dynProbX <- dNmixtureAD_BBNB_s(x, lambda, theta, prob, s, Nmin = Nmin, Nmax = Nmax, len)
              dynCorrectProbX <- 0
              for (N in Nmin:Nmax) {
                dynCorrectProbX <- dynCorrectProbX + dnbinom(N, size = r, prob = pNB) *
                  prod(dBetaBinom_s(x, N, shape1 = alpha, shape2 = beta))
              }
              expect_equal(dynProbX, dynCorrectProbX)
            }
            Nmin <- 0
            Nmax <- 250

            CprobX <- CdNmixtureAD_BBNB_s(x, lambda, theta, prob, s, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixtureAD_BBNB_s(x, lambda, theta, prob, s, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            m <- nimbleModel(code = nc,
                             data = list(x = x),
                             inits = list(lambda = lambda,
                                          prob = prob,
                                          s = s, theta = theta),
                             constants = list(Nmin = Nmin, Nmax = Nmax,
                                              len = len))
            m$calculate()
            MlProbX <- m$getLogProb("x")
            expect_equal(MlProbX, lProbX)

            # Compiled model
            cm <- compileNimble(m)
            cm$calculate()
            CMlProbX <- cm$getLogProb("x")
            expect_equal(CMlProbX, lProbX)

            # All missing
            x <- as.numeric(rep(NA,5))
            probX <- dNmixtureAD_BBNB_s(x, lambda, theta, prob, s, Nmin, Nmax, len)
            expect_equal(probX, 1)

            CprobX <- CdNmixtureAD_BBNB_s(x, lambda, theta, prob, s, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            # Test imputing value for all NAs
            xNA <- c(NA, NA, NA, NA, NA)
            mNA <- nimbleModel(nc, data = list(x = xNA),
                               inits = list(lambda = lambda,
                                            s = s, theta = theta,
                                            prob = prob),
                               constants = list(Nmin = Nmin, Nmax = Nmax,
                                                len = len))

            mNAConf <- configureMCMC(mNA)
            mNAConf$addMonitors('x')
            mNA_MCMC <- buildMCMC(mNAConf)
            cmNA <- compileNimble(mNA, mNA_MCMC)

            set.seed(0)
            cmNA$mNA_MCMC$run(10)

            # Did the imputed values come back?
            expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x[2]"])))

            # Test simulation code
            nSim <- 10
            xSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for (i in 1:nSim) {
              xSim[i,] <- rNmixtureAD_BBNB_s(1, lambda, theta, prob, s, Nmin, Nmax, len)
            }

            CrNmixtureAD_BBNB_s <- compileNimble(rNmixtureAD_BBNB_s)
            CxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i,] <- CrNmixtureAD_BBNB_s(1, lambda, theta, prob, s, Nmin, Nmax, len)
            }
            expect_equal(xSim, CxSim)

            simNodes <- m$getDependencies(c('prob', 'lambda'), self = FALSE)
            mxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for(i in 1:nSim) {
              m$simulate(simNodes, includeData = TRUE)
              mxSim[i,] <- m$x
            }
            expect_equal(mxSim, xSim)

            CmxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for(i in 1:nSim) {
              cm$simulate(simNodes, includeData = TRUE)
              CmxSim[i,] <- cm$x
            }
            expect_identical(CmxSim, mxSim)
          })


# -----------------------------------------------------------------------------
#### 11. Test dNmixtureAD_BBNB_oneObs ####
test_that("dNmixtureAD_BBNB_oneObs works",
          {
            # Uncompiled calculation
            x <- c(1)
            lambda <- 8
            theta <- 1.5
            s <- 2
            prob <- c(0.5)
            Nmin <- max(x)
            Nmax <- 250
            len <- 1

            probX <- dNmixtureAD_BBNB_oneObs(x, lambda, theta, prob, s, Nmin, Nmax)

            # Manually calculate the correct answer
            alpha <- prob * s
            beta <- s - prob * s
            r <- 1 / theta
            pNB <- 1 / (1 + theta * lambda)
            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dBetaBinom_s(x, N, alpha, beta))
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixtureAD_BBNB_oneObs(x, lambda, theta, prob, s, Nmin, Nmax, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Other Nmin / Nmax
            Nmin <- 3
            for(Nmax in 3:6)  {
              dynProbX <- dNmixtureAD_BBNB_oneObs(x, lambda, theta, prob, s, Nmin = Nmin, Nmax = Nmax)
              dynCorrectProbX <- 0
              for (N in Nmin:Nmax) {
                dynCorrectProbX <- dynCorrectProbX + dnbinom(N, size = r, prob = pNB) *
                   prod(dBetaBinom_s(x, N, alpha, beta))
              }
              expect_equal(dynProbX, dynCorrectProbX)
            }
            Nmin <- 0
            Nmax <- 250

            # Compilation and compiled calculations
            call_dNmixtureAD_BBNB_oneObs <- nimbleFunction(
              name = "t11",
              run=function(x=double(),
                           lambda=double(),
                           theta=double(),
                           prob=double(),
                           s=double(),
                           Nmin = double(0),
                           Nmax = double(0),
                           log = integer(0, default=0)) {
                return(dNmixtureAD_BBNB_oneObs(x,lambda,theta,prob,s,Nmin,Nmax,log))
                returnType(double())
              }
            )# Compilation and compiled calculations
            CdNmixtureAD_BBNB_oneObs <- compileNimble(call_dNmixtureAD_BBNB_oneObs)
            CprobX <- CdNmixtureAD_BBNB_oneObs(x, lambda, theta, prob, s, Nmin, Nmax)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixtureAD_BBNB_oneObs(x, lambda, theta, prob, s, Nmin, Nmax, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Dynamic Nmin / Nmax isn't allowed
            expect_error({
              dNmixtureAD_BBNB_oneObs(x, lambda, theta, prob, s, Nmin = -1, Nmax = -1)
            })
            expect_error({
              dNmixtureAD_BBNB_oneObs(x, lambda, theta, prob, s, Nmin, Nmax = -1)
            })
            expect_error({
              dNmixtureAD_BBNB_oneObs(x, lambda, theta, prob, s, Nmin = -1, Nmax)
            })
            expect_error({
              CdNmixtureAD_BBNB_oneObs(x, lambda, theta, prob, s, Nmin = -1, Nmax = -1)
            })
            expect_error({
              CdNmixtureAD_BBNB_oneObs(x, lambda, theta, prob, s, Nmin, Nmax = -1)
            })
            expect_error({
              CdNmixtureAD_BBNB_oneObs(x, lambda, theta, prob, s, Nmin = -1, Nmax)
            })

            # Use in Nimble model
            nc <- nimbleCode({
              x ~ dNmixtureAD_BBNB_oneObs(lambda = lambda, prob = prob,
                                  s = s, theta = theta,
                                  Nmin = Nmin, Nmax = Nmax)
            })

            m <- nimbleModel(code = nc,
                             data = list(x = x),
                             inits = list(lambda = lambda,
                                          prob = prob,
                                          theta = theta,
                                          s = s),
                             constants = list(Nmin = Nmin, Nmax = Nmax))
            m$calculate()
            MlProbX <- m$getLogProb("x")
            expect_equal(MlProbX, lProbX)

            # Compiled model
            cm <- compileNimble(m)
            cm$calculate()
            CMlProbX <- cm$getLogProb("x")
            expect_equal(CMlProbX, lProbX)

            # Missing value
            x <- as.numeric(NA)
            probX <- dNmixtureAD_BBNB_oneObs(x, lambda, theta, prob, s, Nmin, Nmax)
            expect_equal(probX, 1)

            CprobX <- CdNmixtureAD_BBNB_oneObs(x, lambda, theta, prob, s, Nmin, Nmax)
            expect_equal(CprobX, probX)


            # Test imputing value for all NAs
            xNA <- NA
            mNA <- nimbleModel(nc, data = list(x = xNA),
                               inits = list(lambda = lambda,
                                            s = s, theta = theta,
                                            prob = prob),
                               constants = list(Nmin = Nmin, Nmax = Nmax))


            mNAConf <- configureMCMC(mNA)
            mNAConf$addMonitors('x')
            mNA_MCMC <- buildMCMC(mNAConf)
            cmNA <- compileNimble(mNA, mNA_MCMC)

            set.seed(0)
            cmNA$mNA_MCMC$run(10)

            # Did the imputed values come back?
            expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x"])))

            # Test simulation code
            nSim <- 10
            xSim <- numeric(nSim)
            set.seed(1)
            for (i in 1:nSim) {
              xSim[i] <- rNmixtureAD_BBNB_oneObs(1, lambda, theta, prob, s, Nmin, Nmax)
            }

            CrNmixtureAD_BBNB_oneObs <- compileNimble(rNmixtureAD_BBNB_oneObs)
            CxSim <- numeric(nSim)
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i] <- CrNmixtureAD_BBNB_oneObs(1, lambda, theta, prob, s, Nmin, Nmax)
            }
            expect_identical(xSim, CxSim)

            simNodes <- m$getDependencies(c('prob', 'lambda'), self = FALSE)
            mxSim <- numeric(nSim)
            set.seed(1)
            for(i in 1:nSim) {
              m$simulate(simNodes, includeData = TRUE)
              mxSim[i] <- m$x
            }
            expect_identical(mxSim, xSim)

            CmxSim <- numeric(nSim)
            set.seed(1)
            for(i in 1:nSim) {
              cm$simulate(simNodes, includeData = TRUE)
              CmxSim[i] <- cm$x
            }
            expect_identical(CmxSim, mxSim)
          })
