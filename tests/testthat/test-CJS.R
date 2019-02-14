# tests for Cormack-Jolly-Seber distributions

library(testthat)
library(nimble)
context("Testing dCJS-related functions.")

test_that("dCJSss works",
          {
            ## Uncompiled calculation
            x <- c(0, 1, 0, 0)
            probSurvive <- 0.6
            probCapture <- 0.4
            probX <- dCJSss(x, probSurvive, probCapture)
            correctProbX <- probSurvive * (1-probCapture) *
              probSurvive * (probCapture) *
              (probSurvive^2 * (1-probCapture)^2 +
                 probSurvive * (1-probCapture) * (1-probSurvive) +
                 (1-probSurvive))
            expect_equal(probX, correctProbX)

            ## log Prob
            lProbX <- dCJSss(x, probSurvive, probCapture, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            ## Compiles
            CdCJSss <- compileNimble(dCJSss)
            CprobX <- CdCJSss(x, probSurvive, probCapture)
          expect_equal(CprobX, probX)

          ClProbX <- CdCJSss(x, probSurvive, probCapture, log = TRUE)
          expect_equal(ClProbX, lProbX)

          ## r function works
          set.seed(2468)
          oneSim <- rCJSss(1, probSurvive, probCapture, len = 4)
          expect_equal(oneSim, c(0, 1, 0, 0))

          ## Use in model
          nc <- nimbleCode({
            x[1:4] ~ dCJSss(probSurvive, probCapture, len = 0)
            probSurvive ~ dunif(0,1)
            probCapture ~ dunif(0,1)
          })
          m <- nimbleModel(nc, data = list(x = x),
                           inits = list(probSurvive = probSurvive,
                                        probCapture = probCapture))
          m$calculate()
          MlProbX <- m$getLogProb("x")
          expect_equal(MlProbX, lProbX)

          cm <- compileNimble(m, showCompilerOutput = TRUE)
          cm$calculate()
          CMlProbX <- cm$getLogProb("x")
          expect_equal(CMlProbX, lProbX)

          set.seed(2468)
          cm$simulate('x')
          expect_equal(cm$x, c(0, 1, 0, 0))
          })

## Other situations to test:
## imputing values for an NA
## Having the entire capture history be NA (except for first)

test_that("dCJSvv works",
          {
            ## Uncompiled calculation
            x <- c(0,1,0,0)
            probSurvive <- c(0.6, 0.5, 0.4, 0.55)
            probCapture <- c(0.45, 0.5, 0.55, 0.6)
            len <- 4
            probX <- dCJSvv(x, probSurvive, probCapture, len)

            correctProbX <-
              probSurvive[1] * (1 - probCapture[1]) *
              probSurvive[2] * (probCapture[2]) *
              ((probSurvive[3] * (1 - probCapture[3]) * probSurvive[4] * (1 - probCapture[4])) +
               (probSurvive[3] * (1 - probCapture[3]) * (1 - probSurvive[4])) +
                (1 - probSurvive[3]))

            expect_equal(probX, correctProbX)

            ## log Prob
            lProbX <- dCJSvv(x, probSurvive, probCapture, log = TRUE)
            lCorrectProbX <- log(correctProbX)
            expect_equal(lProbX, lCorrectProbX)

            ## Compiles
            CdCJSvv <- compileNimble(dCJSvv)
            CprobX <- CdCJSvv(x, probSurvive, probCapture)
            expect_equal(CprobX, probX)

            ClProbX <- CdCJSvv(x, probSurvive, probCapture, len = 4, log = TRUE)
            expect_equal(ClProbX, lProbX)

            ## r function works
            set.seed(2468)
            oneSim <- rCJSvv(1, probSurvive, probCapture, len = 4)
            expect_equal(oneSim, c(1, 0, 0, 0))

            ## Use in model
            nc <- nimbleCode({
              x[1:4] ~ dCJSvv(probSurvive[1:4], probCapture[1:4], len = 0)
              for (i in 1:4) {
                probSurvive[i] ~ dunif(0,1)
                probCapture[i] ~ dunif(0,1)
              }
            })
            m <- nimbleModel(nc, data = list(x = x),
                             inits = list(probSurvive = probSurvive,
                                          probCapture = probCapture))
            m$calculate()
            MlProbX <- m$getLogProb("x")
            expect_equal(MlProbX, lProbX)

            cm <- compileNimble(m, showCompilerOutput = TRUE)
            cm$calculate()
            CMlProbX <- cm$getLogProb("x")
            expect_equal(CMlProbX, lProbX)

            set.seed(2468)
            cm$simulate('x')
            expect_equal(cm$x, c(0, 1, 0, 0))
          })
