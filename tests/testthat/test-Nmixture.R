# Test the N-mixture distribution nimbleFunction.

# -----------------------------------------------------------------------------
# 0. Load
# Set the context for testthat
context("Testing dNmixture-related functions.")

# -----------------------------------------------------------------------------
# 1. Test dNmixture_v
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
# 2. Test dNmixture_s
test_that("dNmixture_s works",
          {
      # Uncompiled calculation
          x <- c(1, 0, 1, 3, 0)
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
