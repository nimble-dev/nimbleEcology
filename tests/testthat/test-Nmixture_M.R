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
          lambda <- 8
          theta  <- 0.5
          p  <- c(0.2, 0.1, 0.05, 0.2, 0.22)
          J  <- 5
          probX <- dNmixture_MNB_v(x = x, lambda = lambda, p = p, theta = theta, J = J)


      # Calaculate likelihood using truncated infinite sum approach
          K <- 1000
          prob <- c(p, 1 - sum(p))

          min_N <- sum(x)
          max_N <- K


          # Manually calculate the correct answer via trunc inf sum
          N_vec <- min_N:max_N
          l_vec <- numeric(length(N_vec))
          for (i in 1:length(N_vec)) {
            # log prob of N
            lp_N <- dnbinom(x = N_vec[i], size = 1/theta,
                            prob = 1/(1 + theta * lambda), log = TRUE)

            # Log prob of data
            lp_x <- dmultinom(c(x, N_vec[i] - sum(x)), size = N_vec[i], prob = prob,
                              log = TRUE)

            l_vec[i] <- exp(lp_N + lp_x)
          }
          correctProbX <- sum(l_vec)


          expect_equal(probX, correctProbX)

      # Uncompiled log probability

          lProbX        <- dNmixture_MNB_v(x = x, lambda = lambda, p = p, theta = theta, J = J, log = TRUE)
          lCorrectProbX <- log(correctProbX)

          expect_equal(lProbX, lCorrectProbX)

      # Compilation and compiled calculations

          CdNmixture_MNB_v <- compileNimble(dNmixture_MNB_v)
          CprobX <- CdNmixture_MNB_v(x = x, lambda = lambda, p = p, theta = theta, J = J)

          expect_equal(CprobX, probX)

          ClProbX <- CdNmixture_MNB_v(x = x, lambda = lambda, p = p, theta = theta, J = J, log = TRUE)
          expect_equal(ClProbX, lProbX)


      # Use in Nimble model
          nc <- nimbleCode({
            x[1:5] ~ dNmixture_MNB_v(lambda = lambda, p = p[1:5], theta = theta,
                                 J = J)

          })

          m <- nimbleModel(code = nc,
                           data = list(x = x),
                           inits = list(lambda = lambda,
                                        p = p, theta = theta),
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
                 inits = list(lambda = lambda,
                              p = p, theta = theta),
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
            xSim[i,] <- rNmixture_MNB_v(1, lambda = lambda, p = p, theta = theta, J = J)
          }

          CrNmixture_MNB_v <- compileNimble(rNmixture_MNB_v)
          CxSim <- array(NA, dim = c(nSim, J))
          set.seed(1)
          for (i in 1:nSim) {
            CxSim[i,] <- CrNmixture_MNB_v(1, lambda = lambda, p = p, theta = theta, J = J)
          }
          expect_equal(xSim, CxSim)

          simNodes <- m$getDependencies(c('p', 'lambda'), self = FALSE)

          mxSim <- array(NA, dim = c(nSim, J))
          set.seed(1)
          for(i in 1:nSim) {
            m$simulate(simNodes, includeData = TRUE)
            mxSim[i,] <- m$x
          }
          expect_equal(mxSim, xSim)

          CmxSim <- array(NA, dim = c(nSim, J))
          set.seed(1)
          for(i in 1:nSim) {
            cm$simulate(simNodes, includeData = TRUE)
            CmxSim[i,] <- cm$x
          }
          expect_equal(CmxSim, mxSim)
        })

# -----------------------------------------------------------------------------
#### 2. Test dNmixture_MNB_s ####
test_that("dNmixture_MNB_s works",
          {
      # Uncompiled calculation
          x  <- c(1, 0, 1, 3, 2)
          lambda <- 8
          theta  <- 0.5
          p  <- 0.05
          J  <- 5

          probX <- dNmixture_MNB_s(x = x, lambda = lambda, p = p, theta = theta, J = J)


      # Calaculate likelihood using truncated infinite sum approach
          K <- 1000
          prob <- c(rep(p, J), 1 - (p*J))

          min_N <- sum(x)
          max_N <- K

      # Manually calculate the correct answer via trunc inf sum
          N_vec <- min_N:max_N
          l_vec <- numeric(length(N_vec))
          for (i in 1:length(N_vec)) {
            # log prob of N
            lp_N <- dnbinom(x = N_vec[i], size = 1/theta,
                            prob = 1/(1 + theta * lambda), log = TRUE)

            # Log prob of data
            lp_x <- dmultinom(c(x, N_vec[i] - sum(x)), size = N_vec[i], prob = prob,
                              log = TRUE)

            l_vec[i] <- exp(lp_N + lp_x)
          }
          correctProbX <- sum(l_vec)

          expect_equal(correctProbX, probX)

      # Uncompiled log probability
          lProbX <- dNmixture_MNB_s(x = x, lambda = lambda, p = p, theta = theta, J = J, log = TRUE)

      # Compilation and compiled calculations
          CdNmixture_MNB_s <- compileNimble(dNmixture_MNB_s)
          CprobX <- CdNmixture_MNB_s(x = x, lambda = lambda, p = p, theta = theta, J = J)
          expect_equal(CprobX, probX)

          ClProbX <- CdNmixture_MNB_s(x = x, lambda = lambda, p = p, theta = theta, J = J, log = TRUE)
          expect_equal(ClProbX, lProbX)

      # Use in Nimble model
          nc <- nimbleCode({
            x[1:5] ~ dNmixture_MNB_s(lambda = lambda, p = p, theta = theta, J = J)

          })

          m <- nimbleModel(code = nc,
                           data = list(x = x),
                           inits = list(lambda = lambda,
                                        p = p, theta = theta),
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
                 inits = list(lambda = lambda,
                              p = p, theta = theta),
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
            xSim[i,] <- rNmixture_MNB_s(1, lambda = lambda, p = p, theta = theta, J = J)
          }

          CrNmixture_MNB_s <- compileNimble(rNmixture_MNB_s)
          CxSim <- array(NA, dim = c(nSim, J))
          set.seed(1)
          for (i in 1:nSim) {
            CxSim[i,] <- CrNmixture_MNB_s(1, lambda = lambda, p = p, theta = theta, J = J)
          }
          expect_equal(xSim, CxSim)

          simNodes <- m$getDependencies(c('p', 'lambda'), self = FALSE)
          mxSim <- array(NA, dim = c(nSim, J))
          set.seed(1)
          for(i in 1:nSim) {
            m$simulate(simNodes, includeData = TRUE)
            mxSim[i,] <- m$x
          }
          expect_equal(mxSim, xSim)

          CmxSim <- array(NA, dim = c(nSim, J))
          set.seed(1)
          for(i in 1:nSim) {
            cm$simulate(simNodes, includeData = TRUE)
            CmxSim[i,] <- cm$x
          }
          expect_equal(CmxSim, mxSim)
        })


# -----------------------------------------------------------------------------
#### 3. Test dNmixture_MP_v ####
test_that("dNmixture_MP_v works",
          {
            # Uncompiled calculation
            x  <- c(1, 0, 1, 3, 0)
            lambda <- 8
            p  <- c(0.4, 0.2, 0.1, 0.1, 0.05)
            J  <- 5
            probX <- dNmixture_MP_v(x = x, lambda = lambda, p = p, J = J)


            # Calaculate likelihood using truncated infinite sum approach
            K <- 1000
            prob <- c(p, 1 - sum(p))

            min_N <- sum(x)
            max_N <- K


            # Manually calculate the correct answer via trunc inf sum
            N_vec <- min_N:max_N
            l_vec <- numeric(length(N_vec))
            for (i in 1:length(N_vec)) {
              # log prob of N
              lp_N <- dpois(x = N_vec[i], lambda = lambda, log = TRUE)

              # Log prob of data
              lp_x <- dmultinom(c(x, N_vec[i] - sum(x)), size = N_vec[i], prob = prob,
                                log = TRUE)

              l_vec[i] <- exp(lp_N + lp_x)
            }
            correctProbX <- sum(l_vec)

            expect_equal(probX, correctProbX)

            # Uncompiled log probability

            lProbX        <- dNmixture_MP_v(x = x, lambda = lambda, p = p, J = J, log = TRUE)
            lCorrectProbX <- log(correctProbX)

            expect_equal(lProbX, lCorrectProbX)

            # Compilation and compiled calculations

            CdNmixture_MP_v <- compileNimble(dNmixture_MP_v)
            CprobX <- CdNmixture_MP_v(x = x, lambda = lambda, p = p, J = J)

            expect_equal(CprobX, probX)

            ClProbX <- CdNmixture_MP_v(x = x, lambda = lambda, p = p, J = J, log = TRUE)
            expect_equal(ClProbX, lProbX)


            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dNmixture_MP_v(lambda = lambda, p = p[1:5], J = J)

            })

            m <- nimbleModel(code = nc,
                             data = list(x = x),
                             inits = list(lambda = lambda,
                                          p = p),
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
                               inits = list(lambda = lambda,
                                            p = p),
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
              xSim[i,] <- rNmixture_MP_v(1, lambda = lambda, p = p, J = J)
            }

            CrNmixture_MP_v <- compileNimble(rNmixture_MP_v)
            CxSim <- array(NA, dim = c(nSim, J))
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i,] <- CrNmixture_MP_v(1, lambda = lambda, p = p, J = J)
            }
            expect_equal(xSim, CxSim)

            simNodes <- m$getDependencies(c('p', 'lambda'), self = FALSE)

            mxSim <- array(NA, dim = c(nSim, J))
            set.seed(1)
            for(i in 1:nSim) {
              m$simulate(simNodes, includeData = TRUE)
              mxSim[i,] <- m$x
            }
            expect_equal(mxSim, CxSim)

            CmxSim <- array(NA, dim = c(nSim, J))
            set.seed(1)
            for(i in 1:nSim) {
              cm$simulate(simNodes, includeData = TRUE)
              CmxSim[i,] <- cm$x
            }
            expect_equal(CmxSim, mxSim)
          })


# -----------------------------------------------------------------------------
#### 4. Test dNmixture_MP_s ####

test_that("dNmixture_MP_s works",
          {
            # Uncompiled calculation
            x  <- c(1, 0, 1, 3, 0)
            lambda <- 8
            p  <- c(0.1)
            J  <- 5
            probX <- dNmixture_MP_s(x = x, lambda = lambda, p = p, J = J)


            # Calaculate likelihood using truncated infinite sum approach
            K <- 1000
            # prob <- numeric(J + 1)
            # for (i in 1:(J)) {
            #   prob[i] <- pow(1 - p, i - 1) * p
            # }
            # prob[J + 1] <- 1 - sum(prob[1:J])
            prob <- c(rep(p, J), 1- p*J)


            min_N <- sum(x)
            max_N <- K


            # Manually calculate the correct answer via trunc inf sum
            N_vec <- min_N:max_N
            l_vec <- numeric(length(N_vec))
            for (i in 1:length(N_vec)) {
              # log prob of N
              lp_N <- dpois(x = N_vec[i], lambda = lambda, log = TRUE)

              # Log prob of data
              lp_x <- dmultinom(c(x, N_vec[i] - sum(x)), size = N_vec[i], prob = prob,
                                log = TRUE)

              l_vec[i] <- exp(lp_N + lp_x)
            }
            correctProbX <- sum(l_vec)

            expect_equal(probX, correctProbX)

            # Uncompiled log probability

            lProbX        <- dNmixture_MP_s(x = x, lambda = lambda, p = p, J = J, log = TRUE)
            lCorrectProbX <- log(correctProbX)

            expect_equal(lProbX, lCorrectProbX)

            # Compilation and compiled calculations

            CdNmixture_MP_s <- compileNimble(dNmixture_MP_s)
            CprobX <- CdNmixture_MP_s(x = x, lambda = lambda, p = p, J = J)

            expect_equal(CprobX, probX)

            ClProbX <- CdNmixture_MP_s(x = x, lambda = lambda, p = p, J = J, log = TRUE)
            expect_equal(ClProbX, lProbX)


            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dNmixture_MP_s(lambda = lambda, p = p, J = J)

            })

            m <- nimbleModel(code = nc,
                             data = list(x = x),
                             inits = list(lambda = lambda,
                                          p = p),
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
                               inits = list(lambda = lambda,
                                            p = p),
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
              xSim[i,] <- rNmixture_MP_s(1, lambda = lambda, p = p, J = J)
            }

            CrNmixture_MP_s <- compileNimble(rNmixture_MP_s)
            CxSim <- array(NA, dim = c(nSim, J))
            set.seed(1)
            for (i in 1:nSim) {
              CxSim[i,] <- CrNmixture_MP_s(1, lambda = lambda, p = p, J = J)
            }
            expect_equal(xSim, CxSim)

            simNodes <- m$getDependencies(c('p', 'lambda'), self = FALSE)

            mxSim <- array(NA, dim = c(nSim, J))
            set.seed(1)
            for(i in 1:nSim) {
              m$simulate(simNodes, includeData = TRUE)
              mxSim[i,] <- m$x
            }
            expect_equal(mxSim, CxSim)

            CmxSim <- array(NA, dim = c(nSim, J))
            set.seed(1)
            for(i in 1:nSim) {
              cm$simulate(simNodes, includeData = TRUE)
              CmxSim[i,] <- cm$x
            }
            expect_equal(CmxSim, mxSim)
          })

