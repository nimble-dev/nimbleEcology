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
      # Manually calculate the correct answer with special-case Nmin and Nmax
            Nmin <- 3
            for(Nmax in 3:5) {
              probX <- dNmixture_v(x, lambda, prob, Nmin, Nmax, len)
              correctProbX <- 0
              for (N in Nmin:Nmax) {
                correctProbX <- correctProbX + dpois(N, lambda) * prod(dbinom(x, N, prob))
              }
              expect_equal(probX, correctProbX)
            }
            Nmin <- 0
            Nmax <- 250

      # Uncompiled log probability
          lProbX <- dNmixture_v(x, lambda, prob, Nmin, Nmax, len, log = TRUE)
            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dpois(N, lambda) * prod(dbinom(x, N, prob))
            }
            lCorrectProbX <- log(correctProbX)
          expect_equal(lProbX, lCorrectProbX)

      # Dynamic Nmin / Nmax
          dynProbX <- dNmixture_v(x, lambda, prob, Nmin = -1, Nmax = -1, len)
          dNmin <- 0; dNmax <- 21
          dynCorrectProbX <- 0
          for (N in dNmin:dNmax) {
            dynCorrectProbX <- dynCorrectProbX + dpois(N, lambda) * prod(dbinom(x, N, prob))
          }
          expect_equal(dynProbX, dynCorrectProbX)


      # Compilation and compiled calculations
          CdNmixture_v <- compileNimble(dNmixture_v)
            CprobX <- CdNmixture_v(x, lambda, prob, Nmin, Nmax, len)
            probX <- dNmixture_v(x, lambda, prob, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

          ClProbX <- CdNmixture_v(x, lambda, prob, Nmin, Nmax, len, log = TRUE)
          lProbX <- dNmixture_v(x, lambda, prob, Nmin, Nmax, len, log=TRUE)
            expect_equal(ClProbX, lProbX)

          CdynProbX <- CdNmixture_v(x, lambda, prob, Nmin = -1, Nmax = -1, len)
          dynCorrectProbX <- dNmixture_v(x, lambda, prob, Nmin, Nmax, len)
          expect_equal(CdynProbX, dynCorrectProbX)

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

      # Missing value handling
          xna <- c(1, 0, NA, 3, 0)
          correctProbXna <- 0
          for (N in Nmin:Nmax) {
            correctProbXna <- correctProbXna + dpois(N, lambda) * prod(dbinom(xna, N, prob), na.rm=TRUE)
          }
          probXna <- dNmixture_v(xna, lambda, prob, Nmin, Nmax, len)
          expect_equal(probXna, correctProbXna)
          CprobXna <- CdNmixture_v(xna, lambda, prob, Nmin, Nmax, len)
          expect_equal(probXna, CprobXna)

          m_na <- nimbleModel(code = nc,
                           data = list(x = xna),
                           inits = list(lambda = lambda,
                                        prob = prob),
                           constants = list(Nmin = Nmin, Nmax = Nmax,
                                            len = len))
          m_na$calculate()
          MlProbX_na <- m_na$getLogProb("x")
          expect_equal(MlProbX_na, log(correctProbXna))

          cm_na <- compileNimble(m_na)
          cm_na$calculate()
          CMlProbX_na <- cm_na$getLogProb("x")
          expect_equal(CMlProbX_na, log(correctProbXna))

          xna <- c(1, NA, NA, NA, NA)
          correctProbXna <- 0
          for (N in Nmin:Nmax) {
            correctProbXna <- correctProbXna + dpois(N, lambda) * prod(dbinom(xna, N, prob), na.rm=TRUE)
          }
          probXna <- dNmixture_v(xna, lambda, prob, Nmin, Nmax, len)
          expect_equal(probXna, correctProbXna)
          CprobXna <- CdNmixture_v(xna, lambda, prob, Nmin, Nmax, len)
          expect_equal(probXna, CprobXna) 

          xna <- as.numeric(c(NA, NA, NA, NA, NA))
          probXna <- dNmixture_v(xna, lambda, prob, Nmin, Nmax, len)
          expect_equal(probXna, 1)
          CprobXna <- CdNmixture_v(xna, lambda, prob, Nmin, Nmax, len)
          expect_equal(probXna, CprobXna)
          expect_equal(CdNmixture_v(xna, lambda, prob, Nmin, Nmax, len, log=TRUE), 0)

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

      # Dynamic Nmin / Nmax
          dynProbX <- dNmixture_s(x, lambda, prob, Nmin = -1, Nmax = -1, len)
          dNmin <- 0; dNmax <- 20
          dynCorrectProbX <- 0
          for (N in dNmin:dNmax) {
            dynCorrectProbX <- dynCorrectProbX + dpois(N, lambda) * prod(dbinom(x, N, prob))
          }
          expect_equal(dynProbX, dynCorrectProbX)

      # Compilation and compiled calculations
          CdNmixture_s <- compileNimble(dNmixture_s)
          CprobX <- CdNmixture_s(x, lambda, prob, Nmin, Nmax, len)
          expect_equal(CprobX, probX)

          ClProbX <- CdNmixture_s(x, lambda, prob, Nmin, Nmax, len, log = TRUE)
          expect_equal(ClProbX, lProbX)

          CdynProbX <- CdNmixture_s(x, lambda, prob, Nmin = -1, Nmax = -1, len)
          expect_equal(CdynProbX, dynCorrectProbX)

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

      # Missing value handling
          xna <- c(1, 0, NA, 3, 2)
          correctProbXna <- 0
          for (N in Nmin:Nmax) {
            correctProbXna <- correctProbXna + dpois(N, lambda) * prod(dbinom(xna, N, prob), na.rm=TRUE)
          }
          probXna <- dNmixture_s(xna, lambda, prob, Nmin, Nmax, len)
          expect_equal(probXna, correctProbXna)
          CprobXna <- CdNmixture_s(xna, lambda, prob, Nmin, Nmax, len)
          expect_equal(probXna, CprobXna)

          m_na <- nimbleModel(code = nc,
                           data = list(x = xna),
                           inits = list(lambda = lambda,
                                        prob = prob),
                           constants = list(Nmin = Nmin, Nmax = Nmax,
                                            len = len))
          m_na$calculate()
          MlProbX_na <- m_na$getLogProb("x")
          expect_equal(MlProbX_na, log(correctProbXna))

          cm_na <- compileNimble(m_na)
          cm_na$calculate()
          CMlProbX_na <- cm_na$getLogProb("x")
          expect_equal(CMlProbX_na, log(correctProbXna))
  
          xna <- c(1, NA, NA, NA, NA)
          correctProbXna <- 0
          for (N in Nmin:Nmax) {
            correctProbXna <- correctProbXna + dpois(N, lambda) * prod(dbinom(xna, N, prob), na.rm=TRUE)
          }
          probXna <- dNmixture_s(xna, lambda, prob, Nmin, Nmax, len)
          expect_equal(probXna, correctProbXna)
          CprobXna <- CdNmixture_s(xna, lambda, prob, Nmin, Nmax, len)
          expect_equal(probXna, CprobXna) 

          xna <- as.numeric(c(NA, NA, NA, NA, NA))
          probXna <- dNmixture_s(xna, lambda, prob, Nmin, Nmax, len)
          expect_equal(probXna, 1)
          CprobXna <- CdNmixture_s(xna, lambda, prob, Nmin, Nmax, len)
          expect_equal(probXna, CprobXna)
          expect_equal(CdNmixture_s(xna, lambda, prob, Nmin, Nmax, len, log=TRUE), 0)

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

            # Dynamic Nmin/Nmax
            dynProbX <- dNmixture_BNB_v(x, lambda, theta, prob, Nmin = -1, Nmax = -1, len)
            dNmin <- 0; dNmax <- 59
            dynCorrectProbX <- 0
            for (N in dNmin:dNmax) {
              dynCorrectProbX <- dynCorrectProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dbinom(x, N, prob))
            }
            expect_equal(dynProbX, dynCorrectProbX)

            # Some special-case Nmax values
            r <- 1 / theta
            pNB <- 1 / (1 + theta * lambda)
            Nmin <- 3
            for(Nmax in 3:6) {
              probX <- dNmixture_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len)
              # Manually calculate the correct answer
              correctProbX <- 0
              for (N in Nmin:Nmax) {
                correctProbX <- correctProbX + dnbinom(N, size = r, prob = pNB) *
                  prod(dbinom(x, N, prob))
              }
              expect_equal(probX, correctProbX)
            }
            Nmax <- 0
            Nmax <- 250

            # Compilation and compiled calculations
            CdNmixture_BNB_v <- compileNimble(dNmixture_BNB_v)
            CprobX <- CdNmixture_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len)
            probX <- dNmixture_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixture_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            lProbX <- dNmixture_BNB_v(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            CdynProbX <- CdNmixture_BNB_v(x, lambda, theta, prob, Nmin = -1,
                                          Nmax = -1, len)
            dynProbX <- dNmixture_BNB_v(x, lambda, theta, prob, Nmin = -1,
                                          Nmax = -1, len)
            expect_equal(CdynProbX, dynCorrectProbX)

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

            # Missing values
            xna <- c(1, 0, NA, 3, 0)
            probXna <- dNmixture_BNB_v(xna, lambda, theta, prob, Nmin, Nmax, len)

            # Manually calculate the correct answer
            correctProbXna <- 0
            for (N in Nmin:Nmax) {
              correctProbXna <- correctProbXna + dnbinom(N, size = r, prob = pNB) *
                prod(dbinom(xna, N, prob), na.rm=TRUE)
            }
            expect_equal(probXna, correctProbXna)

            # Check compiled version
            CprobXna <- CdNmixture_BNB_v(xna, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(CprobXna, correctProbXna)

            m_na <- nimbleModel(code = nc,
                             data = list(x = xna),
                             inits = list(lambda = lambda,
                                          prob = prob,
                                          theta = theta),
                             constants = list(Nmin = Nmin, Nmax = Nmax,
                                              len = len))
            m_na$calculate()
            MlProbXna <- m_na$getLogProb("x")
            expect_equal(MlProbXna, log(correctProbXna))

            # Compiled model
            cm_na <- compileNimble(m_na)
            cm_na$calculate()
            CMlProbXna <- cm_na$getLogProb("x")
            expect_equal(CMlProbXna, log(correctProbXna))

            xna <- c(1, NA, NA, NA, NA)
            probXna <- dNmixture_BNB_v(xna, lambda, theta, prob, Nmin, Nmax, len)
            # Manually calculate the correct answer
            correctProbXna <- 0
            for (N in Nmin:Nmax) {
              correctProbXna <- correctProbXna + dnbinom(N, size = r, prob = pNB) *
                prod(dbinom(xna, N, prob), na.rm=TRUE)
            }
            expect_equal(probXna, correctProbXna)
            CprobXna <- CdNmixture_BNB_v(xna, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(CprobXna, correctProbXna)

            xna <- as.numeric(rep(NA, 5))
            probXna <- dNmixture_BNB_v(xna, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(probXna, 1)
            CprobXna <- CdNmixture_BNB_v(xna, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(CprobXna, 1)

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

            # Dynamic Nmin/Nmax
            dynProbX <- dNmixture_BNB_s(x, lambda, theta, prob, Nmin = -1, Nmax = -1, len)
            dNmin <- 0; dNmax <- 59
            dynCorrectProbX <- 0
            for (N in dNmin:dNmax) {
              dynCorrectProbX <- dynCorrectProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dbinom(x, N, prob))
            }
            expect_equal(dynProbX, dynCorrectProbX)

            # Compilation and compiled calculations
            CdNmixture_BNB_s <- compileNimble(dNmixture_BNB_s)
            CprobX <- CdNmixture_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixture_BNB_s(x, lambda, theta, prob, Nmin, Nmax, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            CdynProbX <- CdNmixture_BNB_s(x, lambda, theta, prob, Nmin = -1,
                                          Nmax = -1, len)
            expect_equal(CdynProbX, dynCorrectProbX)

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

            # Test NA handling
            xna <- c(1, 0, NA, 3, 0)
            probXna <- dNmixture_BNB_s(xna, lambda, theta, prob, Nmin, Nmax, len)

            # Manually calculate the correct answer
            correctProbXna <- 0
            for (N in Nmin:Nmax) {
              correctProbXna <- correctProbXna + dnbinom(N, size = r, prob = pNB) *
                prod(dbinom(xna, N, prob), na.rm=TRUE)
            }
            expect_equal(probXna, correctProbXna)

            # Check compiled version
            CprobXna <- CdNmixture_BNB_s(xna, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(CprobXna, correctProbXna)

            m_na <- nimbleModel(code = nc,
                             data = list(x = xna),
                             inits = list(lambda = lambda,
                                          prob = prob,
                                          theta = theta),
                             constants = list(Nmin = Nmin, Nmax = Nmax,
                                              len = len))
            m_na$calculate()
            MlProbXna <- m_na$getLogProb("x")
            expect_equal(MlProbXna, log(correctProbXna))

            # Compiled model
            cm_na <- compileNimble(m_na)
            cm_na$calculate()
            CMlProbXna <- cm_na$getLogProb("x")
            expect_equal(CMlProbXna, log(correctProbXna))

            xna <- c(1, NA, NA, NA, NA)
            probXna <- dNmixture_BNB_s(xna, lambda, theta, prob, Nmin, Nmax, len)
            # Manually calculate the correct answer
            correctProbXna <- 0
            for (N in Nmin:Nmax) {
              correctProbXna <- correctProbXna + dnbinom(N, size = r, prob = pNB) *
                prod(dbinom(xna, N, prob), na.rm=TRUE)
            }
            expect_equal(probXna, correctProbXna)
            CprobXna <- CdNmixture_BNB_s(xna, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(CprobXna, correctProbXna)

            xna <- as.numeric(rep(NA, 5))
            probXna <- dNmixture_BNB_s(xna, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(probXna, 1)
            CprobXna <- CdNmixture_BNB_s(xna, lambda, theta, prob, Nmin, Nmax, len)
            expect_equal(CprobXna, 1)

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

            probX <- dNmixture_BNB_oneObs(x, lambda, theta, prob, Nmin, Nmax)

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
            lProbX <- dNmixture_BNB_oneObs(x, lambda, theta, prob, Nmin, Nmax, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Dynamic Nmin/Nmax
            dynProbX <- dNmixture_BNB_oneObs(x, lambda, theta, prob, Nmin = -1, Nmax = -1)
            dNmin <- 1; dNmax <- 17
            dynCorrectProbX <- 0
            for (N in dNmin:dNmax) {
              dynCorrectProbX <- dynCorrectProbX + dnbinom(N, size = r, prob = pNB) *
                prod(dbinom(x, N, prob))
            }
            expect_equal(dynProbX, dynCorrectProbX)

            CdynProbX <- dNmixture_BNB_oneObs(x, lambda, theta, prob, Nmin = -1,
                                              Nmax = -1)
            expect_equal(CdynProbX, dynCorrectProbX)

            # Compilation and compiled calculations
            CdNmixture_BNB_oneObs <- compileNimble(dNmixture_BNB_oneObs)
            CprobX <- CdNmixture_BNB_oneObs(x, lambda, theta, prob, Nmin, Nmax)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixture_BNB_oneObs(x, lambda, theta, prob, Nmin, Nmax, log = TRUE)
            expect_equal(ClProbX, lProbX)

            CdynProbX <- CdNmixture_BNB_oneObs(x, lambda, theta, prob, Nmin = -1, Nmax = -1)
            expect_equal(CdynProbX, dynCorrectProbX)


            # Use in Nimble model
            nc <- nimbleCode({
              x ~ dNmixture_BNB_oneObs(lambda = lambda, prob = prob,
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

            # Test NA handling
            xna <- as.numeric(c(NA))
            probXna <- dNmixture_BNB_oneObs(xna, lambda, theta, prob, Nmin, Nmax)

            # Manually calculate the correct answer
            expect_equal(probXna, 1)

            # Check compiled version
            CprobXna <- CdNmixture_BNB_oneObs(xna, lambda, theta, prob, Nmin, Nmax)
            expect_equal(CprobXna, 1)

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
              xSim[i] <- rNmixture_BNB_oneObs(1, lambda, theta, prob, Nmin, Nmax)
            }

            CrNmixture_BNB_oneObs <- compileNimble(rNmixture_BNB_oneObs)
            CxSim <- numeric(nSim)
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i] <- CrNmixture_BNB_oneObs(1, lambda, theta, prob, Nmin, Nmax)
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
                prod(dBetaBinom_v(x, N, shape1 = alpha, shape2 = beta))
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

            # Dynamic Nmin / Nmax isn't allowed
            expect_error({
              dNmixture_BBP_v(x, lambda, prob, s, Nmin = -1, Nmax = -1, len)
            })
            expect_error({
              dNmixture_BBP_v(x, lambda, prob, s, Nmin = -1, Nmax, len)
            })
            expect_error({
              dNmixture_BBP_v(x, lambda, prob, s, Nmin, Nmax = -1, len)
            })
            expect_error({
              CdNmixture_BBP_v(x, lambda, prob, s, Nmin = -1, Nmax = -1, len)
            })
            expect_error({
              CdNmixture_BBP_v(x, lambda, prob, s, Nmin = -1, Nmax, len)
            })
            expect_error({
              CdNmixture_BBP_v(x, lambda, prob, s, Nmin, Nmax = -1, len)
            })


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

            # Missing values
            # Uncompiled calculation
            x <- c(1, 0, NA, 3, 0)
            probX <- dNmixture_BBP_v(x, lambda, prob, s, Nmin, Nmax, len)

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
            lProbX <- dNmixture_BBP_v(x, lambda, prob, s, Nmin, Nmax, len, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Compilation and compiled calculations
            CdNmixture_BBP_v <- compileNimble(dNmixture_BBP_v)
            CprobX <- CdNmixture_BBP_v(x, lambda, prob, s, Nmin, Nmax, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixture_BBP_v(x, lambda, prob, s, Nmin, Nmax, len, log = TRUE)
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
            probX <- dNmixture_BBP_s(x, lambda, prob, s, Nmin, Nmax, len)
            expect_equal(probX, 1)

            # Compilation and compiled calculations
            CprobX <- CdNmixture_BBP_s(x, lambda, prob, s, Nmin, Nmax, len)
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
                prod(dBetaBinom_s(x, N,
                                  alpha, shape2 = beta))
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


            # Dynamic Nmin / Nmax isn't allowed
            expect_error({
              dNmixture_BBP_s(x, lambda, prob, s, Nmin = -1, Nmax = -1, len)
            })
            expect_error({
              dNmixture_BBP_s(x, lambda, prob, s, Nmin = -1, Nmax, len)
            })
            expect_error({
              dNmixture_BBP_s(x, lambda, prob, s, Nmin, Nmax = -1, len)
            })
            expect_error({
              CdNmixture_BBP_s(x, lambda, prob, s, Nmin = -1, Nmax = -1, len)
            })
            expect_error({
              CdNmixture_BBP_s(x, lambda, prob, s, Nmin = -1, Nmax, len)
            })
            expect_error({
              CdNmixture_BBP_s(x, lambda, prob, s, Nmin, Nmax = -1, len)
            })


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

            # Missing values
            # Uncompiled calculation
            x <- c(1, 0, NA, 3, 0)
            probX <- dNmixture_BBP_s(x, lambda, prob, s, Nmin, Nmax, len)

            # Manually calculate the correct answer
            alpha <- prob * s
            beta <- s - prob * s

            correctProbX <- 0
            for (N in Nmin:Nmax) {
              correctProbX <- correctProbX + dpois(N, lambda) *
                prod(dBetaBinom_s(x, N,
                                  alpha, shape2 = beta), na.rm=TRUE)
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
            probX <- dNmixture_BBP_s(x, lambda, prob, s, Nmin, Nmax, len)
            expect_equal(probX, 1)

            # Compilation and compiled calculations
            CprobX <- CdNmixture_BBP_s(x, lambda, prob, s, Nmin, Nmax, len)
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

            probX <- dNmixture_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax)

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
            lProbX <- dNmixture_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Compilation and compiled calculations
            CdNmixture_BBP_oneObs <- compileNimble(dNmixture_BBP_oneObs)
            CprobX <- CdNmixture_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixture_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Dynamic Nmin / Nmax isn't allowed
            expect_error({
              dNmixture_BBP_oneObs(x, lambda, prob, s, Nmin = -1, Nmax = -1)
            })
            expect_error({
              dNmixture_BBP_oneObs(x, lambda, prob, s, Nmin = -1, Nmax)
            })
            expect_error({
              dNmixture_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax = -1)
            })
            expect_error({
              CdNmixture_BBP_oneObs(x, lambda, prob, s, Nmin = -1, Nmax = -1)
            })
            expect_error({
              CdNmixture_BBP_oneObs(x, lambda, prob, s, Nmin = -1, Nmax)
            })
            expect_error({
              CdNmixture_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax = -1)
            })

            # Use in Nimble model
            nc <- nimbleCode({
              x ~ dNmixture_BBP_oneObs(lambda = lambda, prob = prob,
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
            x <- as.numeric(NA)
            probX <- dNmixture_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax)
            expect_equal(probX, 1)

            # Compilation and compiled calculations
            CprobX <- CdNmixture_BBP_oneObs(x, lambda, prob, s, Nmin, Nmax)
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
              xSim[i] <- rNmixture_BBP_oneObs(1, lambda, prob, s, Nmin, Nmax)
            }

            CrNmixture_BBP_oneObs <- compileNimble(rNmixture_BBP_oneObs)
            CxSim <- numeric(nSim)
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i] <- CrNmixture_BBP_oneObs(1, lambda, prob, s, Nmin, Nmax)
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
                prod(dBetaBinom_v(x, N, shape1 = alpha, shape2 = beta))
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

            # Dynamic Nmin / Nmax isn't allowed
            expect_error({
              dNmixture_BBNB_v(x, lambda, theta, prob, s, Nmin = -1, Nmax = -1, len)
            })
            expect_error({
              dNmixture_BBNB_v(x, lambda, theta, prob, s, Nmin, Nmax = -1, len)
            })
            expect_error({
              dNmixture_BBNB_v(x, lambda, theta, prob, s, Nmin = -1, Nmax, len)
            })
            expect_error({
              CdNmixture_BBNB_v(x, lambda, theta, prob, s, Nmin = -1, Nmax = -1, len)
            })
            expect_error({
              CdNmixture_BBNB_v(x, lambda, theta, prob, s, Nmin, Nmax = -1, len)
            })
            expect_error({
              CdNmixture_BBNB_v(x, lambda, theta, prob, s, Nmin = -1, Nmax, len)
            })

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
                prod(dBetaBinom_s(x, N, shape1 = alpha, shape2 = beta))
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

            # Dynamic Nmin / Nmax isn't allowed
            expect_error({
              dNmixture_BBNB_s(x, lambda, theta, prob, s, Nmin = -1, Nmax = -1, len)
            })
            expect_error({
              dNmixture_BBNB_s(x, lambda, theta, prob, s, Nmin, Nmax = -1, len)
            })
            expect_error({
              dNmixture_BBNB_s(x, lambda, theta, prob, s, Nmin = -1, Nmax, len)
            })
            expect_error({
              CdNmixture_BBNB_s(x, lambda, theta, prob, s, Nmin = -1, Nmax = -1, len)
            })
            expect_error({
              CdNmixture_BBNB_s(x, lambda, theta, prob, s, Nmin, Nmax = -1, len)
            })
            expect_error({
              CdNmixture_BBNB_s(x, lambda, theta, prob, s, Nmin = -1, Nmax, len)
            })
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
            len <- 1

            probX <- dNmixture_BBNB_oneObs(x, lambda, theta, prob, s, Nmin, Nmax)

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
            lProbX <- dNmixture_BBNB_oneObs(x, lambda, theta, prob, s, Nmin, Nmax, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Compilation and compiled calculations
            CdNmixture_BBNB_oneObs <- compileNimble(dNmixture_BBNB_oneObs)
            CprobX <- CdNmixture_BBNB_oneObs(x, lambda, theta, prob, s, Nmin, Nmax)
            expect_equal(CprobX, probX)

            ClProbX <- CdNmixture_BBNB_oneObs(x, lambda, theta, prob, s, Nmin, Nmax, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Dynamic Nmin / Nmax isn't allowed
            expect_error({
              dNmixture_BBNB_oneObs(x, lambda, theta, prob, s, Nmin = -1, Nmax = -1)
            })
            expect_error({
              dNmixture_BBNB_oneObs(x, lambda, theta, prob, s, Nmin, Nmax = -1)
            })
            expect_error({
              dNmixture_BBNB_oneObs(x, lambda, theta, prob, s, Nmin = -1, Nmax)
            })
            expect_error({
              CdNmixture_BBNB_oneObs(x, lambda, theta, prob, s, Nmin = -1, Nmax = -1)
            })
            expect_error({
              CdNmixture_BBNB_oneObs(x, lambda, theta, prob, s, Nmin, Nmax = -1)
            })
            expect_error({
              CdNmixture_BBNB_oneObs(x, lambda, theta, prob, s, Nmin = -1, Nmax)
            })

            # Use in Nimble model
            nc <- nimbleCode({
              x ~ dNmixture_BBNB_oneObs(lambda = lambda, prob = prob,
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
              xSim[i] <- rNmixture_BBNB_oneObs(1, lambda, theta, prob, s, Nmin, Nmax)
            }

            CrNmixture_BBNB_oneObs <- compileNimble(rNmixture_BBNB_oneObs)
            CxSim <- numeric(nSim)
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i] <- CrNmixture_BBNB_oneObs(1, lambda, theta, prob, s, Nmin, Nmax)
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
