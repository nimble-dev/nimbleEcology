

test_that("dBetaBinom works",
          {
            # Uncompiled calculation
            x <- c(4, 2, 8, 0, 3)
            N <- 10
            shape1 <- c(0.1, 2.2, 1.4, 0.4, 0.9)
            shape2 <- c(0.3, 0.5, 0.1, 4, 1.2)
            probX <- dBetaBinom(x, N, shape1, shape2)

            # Manually calculate the correct answer
            correctProbX <- prod(
              choose(N, x) * beta(x + shape1, N - x + shape2) /
                              beta(shape1, shape2)
            )

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dBetaBinom(x, N, shape1, shape2, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Compilation and compiled calculations
            CdBetaBinom <- compileNimble(dBetaBinom)
            CprobX <- CdBetaBinom(x, N, shape1, shape2)
            expect_equal(CprobX, probX)

            ClProbX <- CdBetaBinom(x, N, shape1, shape2, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dBetaBinom(N, shape1[1:5], shape2[1:5])
              N ~ dpois(10)

              for (i in 1:5) {
                shape1[i] ~ dunif(0,10)
                shape2[i] ~ dunif(0,10)
              }
            })
            m <- nimbleModel(nc, data = list(x = x),
                             inits = list(N = N, shape1 = shape1, shape2 = shape2))
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
                               inits = list(N = N, shape1 = shape1, shape2 = shape2))
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
              xSim[i,] <- rBetaBinom(1, N = N, shape1 = shape1, shape2 = shape2)
            set.seed(1)
            CrBetaBinom <- compileNimble(rBetaBinom)
            CxSim <- array(NA, dim = c(nSim, length(x)))
            for(i in 1:nSim)
              CxSim[i,] <- CrBetaBinom(1, N = N, shape1 = shape1, shape2 = shape2)
            expect_equal(xSim, CxSim)

            simNodes <- m$getDependencies(c('shape1', 'shape2', 'N'), self = FALSE)
            mxSim <- array(NA, dim = c(nSim, length(x)))
            set.seed(1)
            for(i in 1:nSim) {
              m$simulate(simNodes, includeData = TRUE)
              mxSim[i,] <- m$x
            }
            expect_equal(mxSim, xSim)

            CmxSim <- array(NA, dim = c(nSim, length(x)))
            set.seed(1)
            for(i in 1:nSim) {
              cm$simulate(simNodes, includeData = TRUE)
              CmxSim[i,] <- cm$x
            }
            expect_equal(CmxSim, mxSim)
          })



test_that("dBetaBinom_One works",
          {
            # Uncompiled calculation
            x <- c(4)
            N <- 10
            shape1 <- c(0.1)
            shape2 <- c(0.3)
            probX <- dBetaBinom_One(x, N, shape1, shape2)

            # Manually calculate the correct answer
            correctProbX <- prod(
              choose(N, x) * beta(x + shape1, N - x + shape2) /
                beta(shape1, shape2)
            )

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dBetaBinom_One(x, N, shape1, shape2, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Compilation and compiled calculations
            CdBetaBinom_One <- compileNimble(dBetaBinom_One)
            CprobX <- CdBetaBinom_One(x, N, shape1, shape2)
            expect_equal(CprobX, probX)

            ClProbX <- CdBetaBinom_One(x, N, shape1, shape2, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Use in Nimble model
            nc <- nimbleCode({
              x ~ dBetaBinom_One(N, shape1, shape2)
              N ~ dpois(10)
              shape1 ~ dunif(0,10)
              shape2 ~ dunif(0,10)
            })
            m <- nimbleModel(nc, data = list(x = x),
                             inits = list(N = N, shape1 = shape1, shape2 = shape2))
            m$calculate()
            MlProbX <- m$getLogProb("x")
            expect_equal(MlProbX, lProbX)

            # Compiled model
            cm <- compileNimble(m)
            cm$calculate()
            CMlProbX <- cm$getLogProb("x")
            expect_equal(CMlProbX, lProbX)

            # Test imputing value for all NAs
            xNA <- c(NA)
            mNA <- nimbleModel(nc, data = list(x = xNA),
                               inits = list(N = N, shape1 = shape1, shape2 = shape2))
            mNAConf <- configureMCMC(mNA)
            mNAConf$addMonitors('x')
            mNA_MCMC <- buildMCMC(mNAConf)
            cmNA <- compileNimble(mNA, mNA_MCMC)

            set.seed(0)
            cmNA$mNA_MCMC$run(10)

            # Did the imputed values come back?
            expect_true(all(!is.na(as.matrix(cmNA$mNA_MCMC$mvSamples)[,"x"])))

            # Test simulation code
            set.seed(1)
            nSim <- 10
            xSim <- array(NA, dim = c(nSim, length(x)))
            for(i in 1:nSim)
              xSim[i,] <- rBetaBinom_One(1, N = N, shape1 = shape1, shape2 = shape2)
            set.seed(1)
            CrBetaBinom_One <- compileNimble(rBetaBinom_One)
            CxSim <- array(NA, dim = c(nSim, length(x)))
            for(i in 1:nSim)
              CxSim[i,] <- CrBetaBinom_One(1, N = N, shape1 = shape1, shape2 = shape2)
            expect_equal(xSim, CxSim)

            simNodes <- m$getDependencies(c('shape1', 'shape2', 'N'), self = FALSE)
            mxSim <- array(NA, dim = c(nSim, length(x)))
            set.seed(1)
            for(i in 1:nSim) {
              m$simulate(simNodes, includeData = TRUE)
              mxSim[i,] <- m$x
            }
            expect_equal(mxSim, xSim)

            CmxSim <- array(NA, dim = c(nSim, length(x)))
            set.seed(1)
            for(i in 1:nSim) {
              cm$simulate(simNodes, includeData = TRUE)
              CmxSim[i,] <- cm$x
            }
            expect_equal(CmxSim, mxSim)
          })
