# Test the N-mixture distribution nimbleFunction.

# -----------------------------------------------------------------------------
# 0. Load
# Set the context for testthat
context("Testing dNmixture-related functions.")

# -----------------------------------------------------------------------------
#### 1. Test dNmixture_v ####
test_that("dNmixture_v works",
          {
      # Uncompiled calculation
          x <- c(1, 0, 1, 3, 0)
          lambda <- 8
          prob <- c(0.5, 0.3, 0.5, 0.4, 0.1)
          Nmin <- 0
          Nmax <- 250
          len <- 5

          probX <- dNmixture_v(x, lambda, prob, Nmin, Nmax, len)
      # Manually calculate the correct answer
          correctProbX <- 0
          for (N in Nmin:Nmax) {
            correctProbX <- correctProbX + dpois(N, lambda) * prod(dbinom(x, N, prob))
          }

          expect_equal(probX, correctProbX)

      # Uncompiled log probability
          lProbX <- dNmixture_v(x, lambda, prob, Nmin, Nmax, len, log = TRUE)
          lCorrectProbX <- log(correctProbX)
          expect_equal(lProbX, lCorrectProbX)

      # Compilation and compiled calculations
          CdNmixture_v <- compileNimble(dNmixture_v)
          CprobX <- CdNmixture_v(x, lambda, prob, Nmin, Nmax, len)
          expect_equal(CprobX, probX)

          ClProbX <- CdNmixture_v(x, lambda, prob, Nmin, Nmax, len, log = TRUE)
          expect_equal(ClProbX, lProbX)

      # Use in Nimble model
          nc <- nimbleCode({
            x[1:5] ~ dNmixture_v(lambda = lambda, prob = prob[1:5],
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
            xSim[i,] <- rNmixture_v(1, lambda, prob, Nmin, Nmax, len)
          }

          CrNmixture_v <- compileNimble(rNmixture_v)
          CxSim <- array(NA, dim = c(nSim, len))
          set.seed(1)
          for (i in 1:nSim) {
            CxSim[i,] <- CrNmixture_v(1, lambda, prob, Nmin, Nmax, len)
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
#### 2. Test dNmixture_s ####
test_that("dNmixture_s works",
          {
      # Uncompiled calculation
          x <- c(1, 0, 1, 3, 2)
          lambda <- 8
          prob <- 0.4
          Nmin <- 0
          Nmax <- 250
          len <- 5

          probX <- dNmixture_s(x, lambda, prob, Nmin, Nmax, len)
      # Manually calculate the correct answer
          correctProbX <- 0
          for (N in Nmin:Nmax) {
            correctProbX <- correctProbX + dpois(N, lambda) * prod(dbinom(x, N, prob))
          }

          expect_equal(probX, correctProbX)

      # Uncompiled log probability
          lProbX <- dNmixture_s(x, lambda, prob, Nmin, Nmax, len, log = TRUE)
          lCorrectProbX <- log(correctProbX)
          expect_equal(lProbX, lCorrectProbX)

      # Compilation and compiled calculations
          CdNmixture_s <- compileNimble(dNmixture_s)
          CprobX <- CdNmixture_s(x, lambda, prob, Nmin, Nmax, len)
          expect_equal(CprobX, probX)

          ClProbX <- CdNmixture_s(x, lambda, prob, Nmin, Nmax, len, log = TRUE)
          expect_equal(ClProbX, lProbX)

      # Use in Nimble model
          nc <- nimbleCode({
            x[1:5] ~ dNmixture_s(lambda = lambda, prob = prob,
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
            xSim[i,] <- rNmixture_s(1, lambda, prob, Nmin, Nmax, len)
          }

          CrNmixture_s <- compileNimble(rNmixture_s)
          CxSim <- array(NA, dim = c(nSim, len))
          set.seed(1)
          for (i in 1:nSim) {
            CxSim[i,] <- CrNmixture_s(1, lambda, prob, Nmin, Nmax, len)
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
#### 3. Test dNmixture_BNB_v ####
test_that("dNmixture_BNB_v works",
          {
            # Uncompiled calculation
            x <- c(1, 0, 3, 3, 0)
            lambda <- 8
            theta <- 2
            prob <- c(0.5, 0.3, 0.5, 0.4, 0.1)
            Nmin <- 0
            Nmax <- 250
            len <- 5

            probX <- dNmixture_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len)

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
            lProbX <- dNmixture_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Compilation and compiled calculations
            CdNmixture_BNB_v <- compileNimble(dNmixture_BNB_v)
            CprobX <- CdNmixture_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixture_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dNmixture_BNB_v(lambda = lambda, prob = prob[1:5],
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
              xSim[i,] <- rNmixture_BNB_v(1, lambda, theta, prob, Nmin, Nmax, len)
            }

            CrNmixture_BNB_v <- compileNimble(rNmixture_BNB_v)
            CxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i,] <- CrNmixture_BNB_v(1, lambda, theta, prob, Nmin, Nmax, len)
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
#### 4. Test dNmixture_BNB_s ####
test_that("dNmixture_BNB_s works",
          {
            # Uncompiled calculation
            x <- c(1, 0, 3, 3, 0)
            lambda <- 8
            theta <- 2
            prob <- 0.4
            Nmin <- 0
            Nmax <- 250
            len <- 5

            probX <- dNmixture_BNB_s(x, lambda, theta = theta, prob, Nmin, Nmax, len)

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
            lProbX <- dNmixture_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Compilation and compiled calculations
            CdNmixture_BNB_s <- compileNimble(dNmixture_BNB_s)
            CprobX <- CdNmixture_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixture_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dNmixture_BNB_s(lambda = lambda, prob = prob,
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
              xSim[i,] <- rNmixture_BNB_s(1, lambda, theta, prob, Nmin, Nmax, len)
            }

            CrNmixture_BNB_s <- compileNimble(rNmixture_BNB_s)
            CxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i,] <- CrNmixture_BNB_s(1, lambda, theta, prob, Nmin, Nmax, len)
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
#### 5. Test dNmixture_BNB_oneObs ####
test_that("dNmixture_BNB_oneObs works",
          {
            # Uncompiled calculation
            x <- c(1)
            lambda <- 8
            theta <- 2
            prob <- c(0.5)
            Nmin <- 0
            Nmax <- 250
            len <- 5

            probX <- dNmixture_BNB_oneObs(x, lambda, theta, prob, Nmin, Nmax, len)

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
            lProbX <- dNmixture_BNB_oneObs(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Compilation and compiled calculations
            CdNmixture_BNB_oneObs <- compileNimble(dNmixture_BNB_oneObs)
            CprobX <- CdNmixture_BNB_oneObs(x, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixture_BNB_oneObs(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Use in Nimble model
            nc <- nimbleCode({
              x ~ dNmixture_BNB_oneObs(lambda = lambda, prob = prob,
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

            # Test imputing value for all NAs
            xNA <- NA
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
            expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x"])))

            # Test simulation code
            nSim <- 10
            xSim <- numeric(nSim)
            set.seed(1)
            for (i in 1:nSim) {
              xSim[i] <- rNmixture_BNB_oneObs(1, lambda, theta, prob, Nmin, Nmax, len)
            }

            CrNmixture_BNB_oneObs <- compileNimble(rNmixture_BNB_oneObs)
            CxSim <- numeric(nSim)
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i] <- CrNmixture_BNB_oneObs(1, lambda, theta, prob, Nmin, Nmax, len)
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
#### 6. Test dNmixture_BBP_v ####
test_that("dNmixture_BBP_v works",
          {
            # Uncompiled calculation
            x <- c(1, 0, 3, 3, 0)
            lambda <- 8
            s <- 2
            prob <- c(0.5, 0.3, 0.5, 0.4, 0.1)
            Nmin <- max(x)
            Nmax <- 250
            len <- 5

            probX <- dNmixture_BBP_v(x, lambda, prob, s, Nmin, Nmax, len)

            # Manually calculate the correct answer
            alpha <- prob * s
            beta <- s - prob * s

            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dpois(N, lambda) *
                prod(dBetaBinom(x, N, alpha = alpha, beta = beta))
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixture_BBP_v(x, lambda, prob, s, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Compilation and compiled calculations
            CdNmixture_BBP_v <- compileNimble(dNmixture_BBP_v)
            CprobX <- CdNmixture_BBP_v(x, lambda, prob, s, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixture_BBP_v(x, lambda, prob, s, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dNmixture_BBP_v(lambda = lambda, prob = prob[1:5],
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
              xSim[i,] <- rNmixture_BBP_v(1, lambda, prob, s, Nmin, Nmax, len)
            }

            CrNmixture_BBP_v <- compileNimble(rNmixture_BBP_v)
            CxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i,] <- CrNmixture_BBP_v(1, lambda, prob, s, Nmin, Nmax, len)
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
#### 7. Test dNmixture_BBP_s ####
test_that("dNmixture_BBP_s works",
          {
            # Uncompiled calculation
            x <- c(1, 0, 3, 3, 0)
            lambda <- 8
            s <- 2
            prob <- 0.4
            Nmin <- max(x)
            Nmax <- 250
            len <- 5

            probX <- dNmixture_BBP_s(x, lambda, prob, s, Nmin, Nmax, len)

            # Manually calculate the correct answer
            alpha <- prob * s
            beta <- s - prob * s

            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dpois(N, lambda) *
                prod(dBetaBinom(x, N,
                                alpha = rep(alpha, len), beta = rep(beta, len)))
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixture_BBP_s(x, lambda, prob, s, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Compilation and compiled calculations
            CdNmixture_BBP_s <- compileNimble(dNmixture_BBP_s)
            CprobX <- CdNmixture_BBP_s(x, lambda, prob, s, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixture_BBP_s(x, lambda, prob, s, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dNmixture_BBP_s(lambda = lambda, prob = prob,
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
              xSim[i,] <- rNmixture_BBP_s(1, lambda, prob, s, Nmin, Nmax, len)
            }

            CrNmixture_BBP_s <- compileNimble(rNmixture_BBP_s)
            CxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i,] <- CrNmixture_BBP_s(1, lambda, prob, s, Nmin, Nmax, len)
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
#### 8. Test dNmixture_BBP_oneObs ####
test_that("dNmixture_BBP_oneObs works",
          {
            # Uncompiled calculation
            x <- c(1)
            lambda <- 8
            s <- 2
            prob <- c(0.5)
            Nmin <- max(x)
            Nmax <- 250
            len <- 5

            probX <- dNmixture_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax, len)

            # Manually calculate the correct answer
            alpha <- prob * s
            beta <- s - prob * s
            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dpois(N, lambda) *
                prod(dBetaBinom_One(x, N, alpha, beta))
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixture_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Compilation and compiled calculations
            CdNmixture_BBP_oneObs <- compileNimble(dNmixture_BBP_oneObs)
            CprobX <- CdNmixture_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixture_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Use in Nimble model
            nc <- nimbleCode({
              x ~ dNmixture_BBP_oneObs(lambda = lambda, prob = prob,
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

            # Test imputing value for all NAs
            xNA <- NA
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
            expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x"])))

            # Test simulation code
            nSim <- 10
            xSim <- numeric(nSim)
            set.seed(1)
            for (i in 1:nSim) {
              xSim[i] <- rNmixture_BBP_oneObs(1, lambda, prob, s, Nmin, Nmax, len)
            }

            CrNmixture_BBP_oneObs <- compileNimble(rNmixture_BBP_oneObs)
            CxSim <- numeric(nSim)
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i] <- CrNmixture_BBP_oneObs(1, lambda, prob, s, Nmin, Nmax, len)
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
#### 9. Test dNmixture_BBNB_v ####
test_that("dNmixture_BBNB_v works",
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

            probX <- dNmixture_BBNB_v(x, lambda, theta, prob, s, Nmin, Nmax, len)

            # Manually calculate the correct answer
            alpha <- prob * s
            beta <- s - prob * s
            r <- 1 / theta
            pNB <- 1 / (1 + theta * lambda)

            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dBetaBinom(x, N, alpha = alpha, beta = beta))
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixture_BBNB_v(x, lambda, theta, prob, s, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Compilation and compiled calculations
            CdNmixture_BBNB_v <- compileNimble(dNmixture_BBNB_v)
            CprobX <- CdNmixture_BBNB_v(x, lambda, theta, prob, s, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixture_BBNB_v(x, lambda, theta, prob, s, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dNmixture_BBNB_v(lambda = lambda, prob = prob[1:5],
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
              xSim[i,] <- rNmixture_BBNB_v(1, lambda, theta, prob, s, Nmin, Nmax, len)
            }

            CrNmixture_BBNB_v <- compileNimble(rNmixture_BBNB_v)
            CxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i,] <- CrNmixture_BBNB_v(1, lambda, theta, prob, s, Nmin, Nmax, len)
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
#### 10. Test dNmixture_BBNB_s ####
test_that("dNmixture_BBNB_s works",
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

            probX <- dNmixture_BBNB_s(x, lambda, theta, prob, s, Nmin, Nmax, len)

            # Manually calculate the correct answer
            alpha <- prob * s
            beta <- s - prob * s
            r <- 1 / theta
            pNB <- 1 / (1 + theta * lambda)

            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dBetaBinom(x, N, alpha = rep(alpha, len), beta = rep(beta, len)))
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixture_BBNB_s(x, lambda, theta, prob, s, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Compilation and compiled calculations
            CdNmixture_BBNB_s <- compileNimble(dNmixture_BBNB_s)
            CprobX <- CdNmixture_BBNB_s(x, lambda, theta, prob, s, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixture_BBNB_s(x, lambda, theta, prob, s, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dNmixture_BBNB_s(lambda = lambda, prob = prob,
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
              xSim[i,] <- rNmixture_BBNB_s(1, lambda, theta, prob, s, Nmin, Nmax, len)
            }

            CrNmixture_BBNB_s <- compileNimble(rNmixture_BBNB_s)
            CxSim <- array(NA, dim = c(nSim, len))
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i,] <- CrNmixture_BBNB_s(1, lambda, theta, prob, s, Nmin, Nmax, len)
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
#### 11. Test dNmixture_BBNB_oneObs ####
test_that("dNmixture_BBNB_oneObs works",
          {
            # Uncompiled calculation
            x <- c(1)
            lambda <- 8
            theta <- 1.5
            s <- 2
            prob <- c(0.5)
            Nmin <- max(x)
            Nmax <- 250
            len <- 5

            probX <- dNmixture_BBNB_oneObs(x, lambda, theta, prob, s, Nmin, Nmax, len)

            # Manually calculate the correct answer
            alpha <- prob * s
            beta <- s - prob * s
            r <- 1 / theta
            pNB <- 1 / (1 + theta * lambda)
            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dBetaBinom_One(x, N, alpha, beta))
            }

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dNmixture_BBNB_oneObs(x, lambda, theta, prob, s, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Compilation and compiled calculations
            CdNmixture_BBNB_oneObs <- compileNimble(dNmixture_BBNB_oneObs)
            CprobX <- CdNmixture_BBNB_oneObs(x, lambda, theta, prob, s, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixture_BBNB_oneObs(x, lambda, theta, prob, s, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Use in Nimble model
            nc <- nimbleCode({
              x ~ dNmixture_BBNB_oneObs(lambda = lambda, prob = prob,
                                  s = s, theta = theta,
                                  Nmin = Nmin, Nmax = Nmax, len = len)

            })

            m <- nimbleModel(code = nc,
                             data = list(x = x),
                             inits = list(lambda = lambda,
                                          prob = prob,
                                          theta = theta,
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

            # Test imputing value for all NAs
            xNA <- NA
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
            expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x"])))

            # Test simulation code
            nSim <- 10
            xSim <- numeric(nSim)
            set.seed(1)
            for (i in 1:nSim) {
              xSim[i] <- rNmixture_BBNB_oneObs(1, lambda, theta, prob, s, Nmin, Nmax, len)
            }

            CrNmixture_BBNB_oneObs <- compileNimble(rNmixture_BBNB_oneObs)
            CxSim <- numeric(nSim)
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i] <- CrNmixture_BBNB_oneObs(1, lambda, theta, prob, s, Nmin, Nmax, len)
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
