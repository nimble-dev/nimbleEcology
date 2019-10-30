# Test the N-mixture distribution nimbleFunction.

# -----------------------------------------------------------------------------
# 0. Load
# Set the context for testthat
context("Testing dNmixture-related functions.")

# -----------------------------------------------------------------------------
# 1. Test dNmixture
test_that("dNmixture works",
          {
      # Uncompiled calculation
          x <- c(1, 0, 1, 3, 0)
          lambda <- 8
          prob <- c(0.5, 0.3, 0.5, 0.4, 0.1)
          minN <- 0
          maxN <- 250
          dynamicMinMax <- 0

          probX <- dNmixture(x, lambda, prob, minN, maxN, dynamicMinMax)
      # Manually calculate the correct answer
          correctProbX <- 0
          for (N in minN:maxN) {
            correctProbX <- correctProbX + dpois(N, lambda) * prod(dbinom(x, N, prob))
          }

          expect_equal(probX, correctProbX)

      # Uncompiled log probability
          lProbX <- dNmixture(x, lambda, prob, minN, maxN, dynamicMinMax, log = TRUE)
          lCorrectProbX <- log(correctProbX)
          expect_equal(lProbX, lCorrectProbX)

      # Compilation and compiled calculations
          CdNmixture <- compileNimble(dNmixture)
          CprobX <- CdNmixture(x, lambda, prob, minN, maxN, dynamicMinMax)
          expect_equal(CprobX, probX)

          ClProbX <- CdNmixture(x, lambda, prob, minN, maxN, dynamicMinMax, log = TRUE)
          expect_equal(ClProbX, lProbX)

      # Use in Nimble model
          nc <- nimbleCode({
            x[1:5] ~ dNmixture(lambda = lambda, prob = prob[1:5],
                               minN = minN, maxN = maxN,
                               dynamicMinMax = dynamicMinMax, len = 5)
          })

          m <- nimbleModel(code = nc,
                           data = list(x = x),
                           inits = list(lambda = lambda,
                                        prob = prob),
                           constants = list(minN = minN, maxN = maxN,
                                            dynamicMinMax = dynamicMinMax))
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
                              prob = prob))
          mNAConf <- configureMCMC(mNA)
          mNAConf$addMonitors('x')
          mNA_MCMC <- buildMCMC(mNAConf)
          cmNA <- compileNimble(mNA, mNA_MCMC)

          set.seed(0)
          cmNA$mNA_MCMC$run(10)

      # Did the imputed values come back?
          expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x[2]"])))

      # Test simulation code
          set.seed(1)
          nSim <- 10
          xSim <- array(NA, dim = c(nSim, length(x)))
          for(i in 1:nSim)
            xSim[i,] <- rCJS_ss(1, probSurvive, probCapture, len = length(x))
          set.seed(1)
          CrCJS_ss <- compileNimble(rCJS_ss)
          CxSim <- array(NA, dim = c(nSim, length(x)))
          for(i in 1:nSim)
            CxSim[i,] <- CrCJS_ss(1, probSurvive, probCapture, len = length(x))
          expect_identical(xSim, CxSim)

          simNodes <- m$getDependencies(c('probSurvive', 'probCapture'), self = FALSE)
          mxSim <- array(NA, dim = c(nSim, length(x)))
          set.seed(1)
          for(i in 1:nSim) {
            m$simulate(simNodes, includeData = TRUE)
            mxSim[i,] <- m$x
          }
          expect_identical(mxSim, xSim)

          CmxSim <- array(NA, dim = c(nSim, length(x)))
          set.seed(1)
          for(i in 1:nSim) {
            cm$simulate(simNodes, includeData = TRUE)
            CmxSim[i,] <- cm$x
          }
          expect_identical(CmxSim, mxSim)
        })
