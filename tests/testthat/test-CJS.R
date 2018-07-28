# tests for Cormack-Jolly-Seber distributions

library(testthat)
context("Testing dCJS-related functions.")

test_that("dCJS_ss works",
          {
            ## Uncompiled calculation
            x <- c(1, 0, 1, 0, 0)
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
          oneSim <- rCJSss(1, probSurvive, probCapture, len = 5)
          expect_equal(oneSim, c(1, 0, 1, 0, 0))

          ## Use in model
          nc <- nimbleCode({
            x[1:5] ~ dCJSss(probSurvive, probCapture, len = 0)
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
          expect_equal(cm$x, c(1, 0, 1, 0, 0))
          })

## Other situations to test:
## imputing values for an NA
## Having the entire capture history be NA (except for first)

test_that("dCJS_vv works",
          {
            ## Uncompiled calculation
            x <- c(1, 0, 1, 0, 0)
            ## STOPPED HERE
            probSurvive <- c(0.6, 0.5, 0.4, 0.55, 0.65)
            probCapture <- c(0.4, 0.45, 0.5, 0.55, 0.6)
            probX <- dCJSss(x, probSurvive, probCapture)
            correctProbX <- probSurvive[1] * (1-probCapture[1]) *
              probSurvive[2] * (probCapture) *
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
            oneSim <- rCJSss(1, probSurvive, probCapture, len = 5)
            expect_equal(oneSim, c(1, 0, 1, 0, 0))

            ## Use in model
            nc <- nimbleCode({
              x[1:5] ~ dCJSss(probSurvive, probCapture, len = 0)
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
            expect_equal(cm$x, c(1, 0, 1, 0, 0))
          })
