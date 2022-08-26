# Test the N-mixture_M distribution nimbleFunction.

# -----------------------------------------------------------------------------
# 0. Load
# Set the context for testthat
context("Testing dNmixture_M-related functions.")

# -----------------------------------------------------------------------------
#### 1. Test dNmixture_MNB_v ####
test_that("dNmixture_MNB_v works",
          {
      # Uncompiled calculation
          x  <- c(1, 0, 1, 3, 0)
          mu <- 8
          r  <- 0.5
          p  <- c(0.5, 0.3, 0.5, 0.4, 0.1)
          J  <- 5
          probX <- dNmixture_MNB_v(x = x, mu = mu, p = p, r = r, J = J)

      # Manually calculate the correct answer

          x_tot <- sum(x)
          x_miss <- sum(x * nimSeq(0, J - 1))
          pp <- nimC(0, p)
          prob <- nimNumeric(J)
          for (j in 1:J) {
              prob[j] <- prod(1 - pp[1:j]) * p[j]
          }
          ptot <- sum(prob)
          term1 <- lgamma(r + x_tot) - lgamma(r) - sum(lfactorial(x))
          term2 <- r * log(r) + x_tot * log(mu)
          term3 <- sum(x * log(prob))
          term4 <- -(x_tot + r) * log(r + mu * ptot)
          logProb <- term1 + term2 + term3 + term4
          correctProbX <- exp(logProb)

          expect_equal(probX, correctProbX)

      # Uncompiled log probability

          lProbX        <- dNmixture_MNB_v(x = x, mu = mu, p = p, r = r, J = J, log = TRUE)
          lCorrectProbX <- log(correctProbX)

          expect_equal(lProbX, lCorrectProbX)

      # Compilation and compiled calculations

          CdNmixture_MNB_v <- compileNimble(dNmixture_MNB_v)
          CprobX <- CdNmixture_MNB_v(x = x, mu = mu, p = p, r = r, J = J)

          expect_equal(CprobX, probX)

          ClProbX <- CdNmixture_MNB_v(x = x, mu = mu, p = p, r = r, J = J, log = TRUE)
          expect_equal(ClProbX, lProbX)


      # Use in Nimble model
          nc <- nimbleCode({
            x[1:5] ~ dNmixture_MNB_v(mu = mu, p = p[1:5], r = r,
                                 J = J)

          })

          m <- nimbleModel(code = nc,
                           data = list(x = x),
                           inits = list(mu = mu,
                                        p = p, r = r),
                           constants = list(J = J))
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
                 inits = list(mu = mu,
                              p = p, r = r),
                 constants = list(J = J))


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
          xSim <- array(NA, dim = c(nSim, J))
          set.seed(1)
          for (i in 1:nSim) {
            xSim[i,] <- rNmixture_MNB_v(1, mu = mu, p = p, r = r, J = J)
          }

          CrNmixture_MNB_v <- compileNimble(rNmixture_MNB_v)
          CxSim <- array(NA, dim = c(nSim, J))
          set.seed(1)
          for (i in 1:nSim) {
            CxSim[i,] <- CrNmixture_MNB_v(1, mu = mu, p = p, r = r, J = J)
          }
          expect_identical(xSim, CxSim)

          simNodes <- m$getDependencies(c('p', 'mu'), self = FALSE)

          mxSim <- array(NA, dim = c(nSim, J))
          set.seed(1)
          for(i in 1:nSim) {
            m$simulate(simNodes, includeData = TRUE)
            mxSim[i,] <- m$x
          }
          expect_identical(mxSim, xSim)

          CmxSim <- array(NA, dim = c(nSim, J))
          set.seed(1)
          for(i in 1:nSim) {
            cm$simulate(simNodes, includeData = TRUE)
            CmxSim[i,] <- cm$x
          }
          expect_identical(CmxSim, mxSim)
        })

# -----------------------------------------------------------------------------
#### 2. Test dNmixture_MNB_s ####
test_that("dNmixture_MNB_s works",
          {
      # Uncompiled calculation
          x  <- c(1, 0, 1, 3, 2)
          mu <- 8
          r  <- 0.5
          p  <- 0.4
          J  <- 5

          probX <- dNmixture_MNB_s(x = x, mu = mu, p = p, r = r, J = J)

      # Uncompiled log probability
          lProbX <- dNmixture_MNB_s(x = x, mu = mu, p = p, r = r, J = J, log = TRUE)

      # Compilation and compiled calculations
          CdNmixture_MNB_s <- compileNimble(dNmixture_MNB_s)
          CprobX <- CdNmixture_MNB_s(x = x, mu = mu, p = p, r = r, J = J)
          expect_equal(CprobX, probX)

          ClProbX <- CdNmixture_MNB_s(x = x, mu = mu, p = p, r = r, J = J, log = TRUE)
          expect_equal(ClProbX, lProbX)

      # Use in Nimble model
          nc <- nimbleCode({
            x[1:5] ~ dNmixture_MNB_s(mu = mu, p = p, r = r, J = J)

          })

          m <- nimbleModel(code = nc,
                           data = list(x = x),
                           inits = list(mu = mu,
                                        p = p, r = r),
                           constants = list(J = J))
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
                 inits = list(mu = mu,
                              p = p, r = r),
                 constants = list(J = J))


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
          xSim <- array(NA, dim = c(nSim, J))
          set.seed(1)
          for (i in 1:nSim) {
            xSim[i,] <- rNmixture_MNB_s(1, mu = mu, p = p, r = r, J = J)
          }

          CrNmixture_MNB_s <- compileNimble(rNmixture_MNB_s)
          CxSim <- array(NA, dim = c(nSim, J))
          set.seed(1)
          for (i in 1:nSim) {
            CxSim[i,] <- CrNmixture_MNB_s(1, mu = mu, p = p, r = r, J = J)
          }
          expect_identical(xSim, CxSim)

          simNodes <- m$getDependencies(c('p', 'mu'), self = FALSE)
          mxSim <- array(NA, dim = c(nSim, J))
          set.seed(1)
          for(i in 1:nSim) {
            m$simulate(simNodes, includeData = TRUE)
            mxSim[i,] <- m$x
          }
          expect_identical(mxSim, xSim)

          CmxSim <- array(NA, dim = c(nSim, J))
          set.seed(1)
          for(i in 1:nSim) {
            cm$simulate(simNodes, includeData = TRUE)
            CmxSim[i,] <- cm$x
          }
          expect_identical(CmxSim, mxSim)
        })
