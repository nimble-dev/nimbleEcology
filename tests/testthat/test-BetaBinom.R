

test_that("dBetaBinom_v works",
          {
            # Uncompiled calculation
            x <- c(4, 2, 8, 0, 3)
            N <- 10
            shape1 <- c(0.1, 2.2, 1.4, 0.4, 0.9)
            shape2 <- c(0.3, 0.5, 0.1, 4, 1.2)
            probX <- dBetaBinom_v(x, N, shape1, shape2)
            len <- 5

            # Manually calculate the correct answer
            correctProbX <- prod(
              choose(N, x) * beta(x + shape1, N - x + shape2) /
                              beta(shape1, shape2)
            )

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dBetaBinom_v(x, N, shape1, shape2, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Compilation and compiled calculations
            call_dBetaBinom_v <- nimbleFunction(
              name = "call_dBetaBinom_v",
              run = function(x = double(1),
                             N = double(),
                             shape1 = double(1),
                             shape2 = double(1),
                             len = double(),
                             log = integer(0, default = 0)) {
                return(dBetaBinom_v(x, N, shape1, shape2, len, log))
                returnType(double())
              }
            )

            # Compilation and compiled calculations
            CdBetaBinom_v <- compileNimble(call_dBetaBinom_v)
            CprobX <- CdBetaBinom_v(x, N, shape1, shape2, len = len)
            expect_equal(CprobX, probX)

            ClProbX <- CdBetaBinom_v(x, N, shape1, shape2, len, log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dBetaBinom_v(N, shape1[1:5], shape2[1:5], len=5)
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

            # Missing values
            xna <- c(4, NA, 8, 0, 3)
            probXna <- dBetaBinom_v(xna, N, shape1, shape2)
            len <- 5

            correctProbXna <- prod(
              choose(N, xna) * beta(xna + shape1, N - xna + shape2) /
                              beta(shape1, shape2)
            , na.rm=TRUE)
            expect_equal(probXna, correctProbXna)

            CprobXna <- CdBetaBinom_v(xna, N, shape1, shape2, len = len)
            expect_equal(CprobXna, probXna)

            # All NAs
            xna <- as.numeric(rep(NA, 5))
            expect_equal(CdBetaBinom_v(xna, N, shape1, shape2, len=len), 1)
            expect_equal(CdBetaBinom_v(xna, N, shape1, shape2, len=len, log=TRUE), 0)

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
              xSim[i,] <- rBetaBinom_v(1, N = N, shape1 = shape1, shape2 = shape2, len=len)
            set.seed(1)
            CrBetaBinom_v <- compileNimble(rBetaBinom_v)
            CxSim <- array(NA, dim = c(nSim, length(x)))
            for(i in 1:nSim)
              CxSim[i,] <- CrBetaBinom_v(1, N = N, shape1 = shape1, shape2 = shape2, len=len)
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



test_that("dBetaBinom_s works",
          {
            # Uncompiled calculation
            x <- c(4, 2, 8, 0, 3)
            N <- 10
            shape1 <- c(0.1)
            shape2 <- c(0.3)
            len <- 5
            probX <- dBetaBinom_s(x, N, shape1, shape2, len)

            # Manually calculate the correct answer
            correctProbX <- prod(
              choose(N, x) * beta(x + shape1, N - x + shape2) /
                beta(shape1, shape2)
            )

            expect_equal(probX, correctProbX)

            # Uncompiled log probability
            lProbX <- dBetaBinom_s(x, N, shape1, shape2, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            # Compilation and compiled calculations
            call_dBetaBinom_s <- nimbleFunction(
                          name = "call_dBetaBinom_s",
                          run = function(x = double(1),
                             N = double(),
                             shape1 = double(),
                             shape2 = double(),
                             len = double(),
                             log = integer(0, default = 0)) {
                return(dBetaBinom_s(x, N, shape1, shape2, len, log))
                returnType(double())
              }
            )

            # Compilation and compiled calculations
            CdBetaBinom_s <- compileNimble(call_dBetaBinom_s)
            CprobX <- CdBetaBinom_s(x, N, shape1, shape2, len)
            expect_equal(CprobX, probX)

            ClProbX <- CdBetaBinom_s(x, N, shape1, shape2, len,log = TRUE)
            expect_equal(ClProbX, lProbX)

            # Use in Nimble model
            nc <- nimbleCode({
              x[1:5] ~ dBetaBinom_s(N, shape1, shape2,len=5)
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

            # Missing values
            xna <- c(4, NA, 8, 0, 3)
            probXna <- dBetaBinom_s(xna, N, shape1, shape2)
            len <- 5

            correctProbXna <- prod(
              choose(N, xna) * beta(xna + shape1, N - xna + shape2) /
                              beta(shape1, shape2)
            , na.rm=TRUE)
            expect_equal(probXna, correctProbXna)

            CprobXna <- CdBetaBinom_s(xna, N, shape1, shape2, len = len)
            expect_equal(CprobXna, probXna)

            # All NAs
            xna <- as.numeric(rep(NA, 5))
            expect_equal(CdBetaBinom_s(xna, N, shape1, shape2, len=len), 1)
            expect_equal(CdBetaBinom_s(xna, N, shape1, shape2, len=len, log=TRUE), 0)

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
              xSim[i,] <- rBetaBinom_s(1, N = N, shape1 = shape1, shape2 = shape2, len)
            set.seed(1)
            CrBetaBinom_s <- compileNimble(rBetaBinom_s)
            CxSim <- array(NA, dim = c(nSim, length(x)))
            for(i in 1:nSim)
              CxSim[i,] <- CrBetaBinom_s(1, N = N, shape1 = shape1, shape2 = shape2, len)
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
