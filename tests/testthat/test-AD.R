# Testing examples:

# install nimble from branch ADoak: devtools::install_github("nimble-dev/nimble", ref = "ADoak", subdir = "packages/nimble")
# install nimbleEcology from branch AD_0.3: devtools::install_github("nimble-dev/nimbleEcology", ref = "AD_0.3")

# load nimble's testing tools
library(nimble)
library(nimbleEcology)
# source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
# source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
source("../nimble/packages/nimble/tests/testthat/test_utils.R")
source("../nimble/packages/nimble/tests/testthat/AD_test_utils.R")

EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)
nimbleOptions(allowDynamicIndexing = FALSE)

#####################
#### dOcc_s case ####

dat <- c(1,1,0,0) # A vector of observations
probOcc <- 0.6
probDetect <- 0.4

nc <- nimbleCode({
  x[1:4] ~ dOcc_s(probOcc, probDetect, len = 4)
  probOcc ~ dunif(0,1)
  probDetect ~ dunif(0,1)
})

Rmodel <- nimbleModel(nc, data = list(x = dat),
                         inits = list(probOcc = probOcc,
                                      probDetect = probDetect),
                      buildDerivs=TRUE)

Cmodel <- compileNimble(Rmodel)

nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c('probOcc', 'probDetect'))
v1_case1 <- list(arg1 = c(0.6, 0.4)) # taping values for probOcc and probDetect
v2_case1 <- list(arg1 = c(0.65, 0.35)) # testing values for probOcc and probDetect
RCrelTol = c(1e-15, 1e-8, 1e-3, 1e-14)

model_calculate_test_case(Rmodel, Cmodel,
                          model_calculate_test, nodesList_case1,
                          v1_case1, v2_case1,
                          0:2) # lots of output numbers with no warning messages means it passes.

#####################
#### dOcc_v case ####

x <- c(1,0,1,1,0)
probOcc <- 0.4
probDetect <- c(0.7, 0.3, 0.5, 0.7, 0.25)

probOcc2 <- 0.5
probDetect2 <- c(0.77, 0.39, 0.52, 0.78, 0.32)


nc <- nimbleCode({
  x[1:5] ~ dOcc_v(probOcc, probDetect[1:5], len = 5)
  for (i in 1:5) {
    probDetect[i] ~ dunif(0,1)
  }
  probOcc ~ dunif(0,1)
})
Rmodel <- nimbleModel(nc, data = list(x = x),
                 inits = list(probOcc = probOcc,
                              probDetect = probDetect),
                 buildDerivs=TRUE)
Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)

nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel,
                                                c('probOcc', Rmodel$expandNodeNames('probDetect[1:5]')))
v1_case1 <- list(arg1 = c(probOcc, probDetect)) # taping values for probOcc and probDetect
v2_case1 <- list(arg1 = c(probOcc2, probDetect2)) # testing values for probOcc and probDetect

model_calculate_test_case(Rmodel, Cmodel,
                          model_calculate_test, nodesList_case1,
                          v1_case1, v2_case1,
                          0:2)


##########################
#### dNmixture_s case ####

x <- c(7, 7, 6, 9, 10)
lambda <- 15
prob <- 0.7

lambda2 <- 18
prob2 <- 0.5


nc <- nimbleCode({
  x[1:5] ~ dNmixture_s(lambda, prob,
                       Nmin = 0, Nmax = 100, len = 5)
  prob ~ dunif(0, 1)
  lambda ~ dunif(0, 100)
})
Rmodel <- nimbleModel(nc, data = list(x = x),
                 inits = list(prob = prob,
                              lambda = lambda),
                 buildDerivs=TRUE)
Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)
Cmodel$calculate()

nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c('prob', 'lambda'))
v1_case1 <- list(arg1 = c(prob, lambda)) # taping values for prob and lambda
v2_case1 <- list(arg1 = c(prob2, lambda2)) # testing values for prob and lambda

model_calculate_test_case(Rmodel, Cmodel,
                          model_calculate_test, nodesList_case1,
                          v1_case1, v2_case1,
                          0:2)


######################
#### dCJS_ss case ####

x <- c(1, 0, 1, 0, 0, 0)
probSurvive <- 0.8
probCapture <- 0.4

probSurvive2 <- 0.7
probCapture2 <- 0.2


nc <- nimbleCode({
  x[1:6] ~ dCJS_ss(probSurvive, probCapture, len = 6)
  probSurvive ~ dunif(0, 1)
  probCapture ~ dunif(0, 1)
})
Rmodel <- nimbleModel(nc, data = list(x = x),
                 inits = list(probSurvive = probSurvive,
                              probCapture = probCapture),
                 buildDerivs=TRUE)
Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)
Cmodel$calculate()

nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c('probSurvive', 'probCapture'))
v1_case1 <- list(arg1 = c(probSurvive, probCapture)) # taping values for prob and lambda
v2_case1 <- list(arg1 = c(probSurvive2, probCapture2)) # testing values for prob and lambda

model_calculate_test_case(Rmodel, Cmodel,
                          model_calculate_test, nodesList_case1,
                          v1_case1, v2_case1,
                          0:2)


######################
#### dCJS_sv case ####

x <- c(1, 0, 1, 0, 0, 0)
probSurvive <- 0.8
probCapture <- c(1, 0.5, 0.5, 0.4, 0.3, 0.4)

probSurvive2 <- 0.7
probCapture2 <- c(1, 0.6, 0.7, 0.4, 0.2, 0.2)


nc <- nimbleCode({
  x[1:6] ~ dCJS_sv(probSurvive, probCapture[1:6], len = 6)
  for (i in 1:6) {
    probCapture[i] ~ dunif(0, 1)
  }
  probSurvive ~ dunif(0, 1)

})
Rmodel <- nimbleModel(nc, data = list(x = x),
                 inits = list(probSurvive = probSurvive,
                              probCapture = probCapture),
                 buildDerivs=TRUE)
Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)
Cmodel$calculate()

nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c('probSurvive', Rmodel$expandNodeNames('probCapture[2:6]')))
v1_case1 <- list(arg1 = c(probSurvive, probCapture[2:6])) # taping values for prob and lambda
v2_case1 <- list(arg1 = c(probSurvive2, probCapture2[2:6])) # testing values for prob and lambda

model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                          nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                          order = 0:2)


######################
#### dCJS_vs case ####

x <- c(1, 0, 1, 0, 0, 0)
probSurvive <- c(0.8, 0.5, 0.3, 0.9, 0.9)
probCapture <- 0.5

probSurvive2 <- c(0.7, 0.55, 0.32, 0.8, 0.1)
probCapture2 <- 0.7


nc <- nimbleCode({
  x[1:6] ~ dCJS_vs(probSurvive[1:5], probCapture, len = 6)
  for (i in 1:5) {
    probSurvive[i] ~ dunif(0, 1)
  }
  probCapture ~ dunif(0, 1)

})
Rmodel <- nimbleModel(nc, data = list(x = x),
                 inits = list(probSurvive = probSurvive,
                              probCapture = probCapture),
                 buildDerivs=TRUE)
Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)
Cmodel$calculate()

nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c(Rmodel$expandNodeNames('probSurvive[1:5]'), 'probCapture'))
v1_case1 <- list(arg1 = c(probSurvive, probCapture)) # taping values for prob and lambda
v2_case1 <- list(arg1 = c(probSurvive2, probCapture2)) # testing values for prob and lambda

model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                          nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                          order = 0:2)

######################
#### dCJS_vv case ####

x <- c(1, 0, 1, 0, 0, 0)
probSurvive <- c(0.8, 0.5, 0.3, 0.9, 0.9)
probCapture <- c(1, 0.5, 0.5, 0.4, 0.3, 0.4)

probSurvive2 <- c(0.7, 0.55, 0.32, 0.8, 0.1)
probCapture2 <- c(1, 0.6, 0.7, 0.4, 0.2, 0.2)


nc <- nimbleCode({
  x[1:6] ~ dCJS_vv(probSurvive[1:5], probCapture[1:6], len = 6)
  for (i in 1:5) {
    probSurvive[i] ~ dunif(0, 1)
  }
  for (i in 1:6) {
    probCapture[i] ~ dunif(0, 1)
  }

})
Rmodel <- nimbleModel(nc, data = list(x = x),
                 inits = list(probSurvive = probSurvive,
                              probCapture = probCapture),
                 buildDerivs=TRUE)
Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)
Cmodel$calculate()

nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c(Rmodel$expandNodeNames('probSurvive[1:5]'),
                                                                       Rmodel$expandNodeNames('probCapture[2:6]')))
v1_case1 <- list(arg1 = c(probSurvive, probCapture[2:6])) # taping values for prob and lambda
v2_case1 <- list(arg1 = c(probSurvive2, probCapture2[2:6])) # testing values for prob and lambda

model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                          nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                          order = 0:2)

######################
#### dHMM case ####

x <- c(1, 1, 1, 2, 2)

init <- c(0.4, 0.2, 0.4)
probObs <- t(array(
         c(0.9, 0.1,
           0.1, 0.9,
           0.8, 0.2),
         c(2, 3)))

probTrans <- t(array(
          c(0.3, 0.4, 0.2,
            0.1, 0.1, 0.8,
            0.05, 0.05, 0.9),
          c(3,3)))

init2 <- c(0.6, 0.1, 0.3)
probObs2 <- t(array(
         c(0.9, 0.1,
           0.1, 0.9,
           0.7, 0.3),
         c(2, 3)))

probTrans2 <- t(array(
          c(0.4, 0.4, 0.2,
            0.05, 0.25, 0.7,
            0.05, 0.15, 0.8),
          c(3,3)))

nc <- nimbleCode({
  x[1:5] ~ dHMM(init[1:3], probObs = probObs[1:3,1:2],
                  probTrans = probTrans[1:3, 1:3], len = 5, checkRowSums = 0)
  for (i in 1:2) {
    init[i] ~ dunif(0, 1)
  }
  init[3] <- 1 - init[1] - init[2]

  for (i in 1:3) {
    probObs[i, 1] ~ dunif(0, 1)
    probObs[i, 2] <- 1 - probObs[i, 1]
    probTrans[i, 1] ~ dunif(0, 1)
    probTrans[i, 2] ~ dunif(0, 1 - probTrans[i, 1])
    probTrans[i, 3] <- 1 - probTrans[i, 1] - probTrans[i, 2]
  }

})
Rmodel <- nimbleModel(nc, data = list(x = x),
                 inits = list(
                   init = init,
                   probObs = probObs,
                   probTrans = probTrans
                 ),
                 buildDerivs=TRUE)
Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)
Cmodel$calculate()

nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c(Rmodel$expandNodeNames('init[1:3]'),
                                                                       Rmodel$expandNodeNames('probObs[1:3, 1:2]'),
                                                                       Rmodel$expandNodeNames('probTrans[1:3, 1:3]')
                                                                       ))
v1_case1 <- list(arg1 = c(init[1:3],  probObs[1:3, 1:2], probTrans[1:3, 1:3]))
v2_case1 <- list(arg1 = c(init2[1:3], probObs2[1:3, 1:2],  probTrans2[1:3, 1:3]))

model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                          nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                          order = 0:2)
#######


######################
#### dHMM with 0s in transition matrix case ####

x <- c(1, 1, 1, 2, 2)

init <- c(0.4, 0.2, 0.4)
probObs <- t(array(
         c(0.9, 0.1,
           0.1, 0.9,
           0.8, 0.2),
         c(2, 3)))

probTrans <- t(array(
          c(0.3, 0.4, 0.2,
            0.1, 0.1, 0.8,
            0.05, 0.05, 0.9),
          c(3,3)))

init2 <- c(0.6, 0.1, 0.3)
probObs2 <- t(array(
         c(1, 0,
           0, 1,
           0.7, 0.3),
         c(2, 3)))

probTrans2 <- t(array(
          c(0.4, 0.4, 0.2,
            0, 0.3, 0.7,
            0, 0, 1),
          c(3,3)))

nc <- nimbleCode({
  x[1:5] ~ dHMM(init[1:3], probObs = probObs[1:3,1:2],
                  probTrans = probTrans[1:3, 1:3], len = 5, checkRowSums = 0)
  for (i in 1:2) {
    init[i] ~ dunif(0, 1)
  }
  init[3] <- 1 - init[1] - init[2]

  for (i in 1:3) {
    probObs[i, 1] ~ dunif(0, 1)
    probObs[i, 2] <- 1 - probObs[i, 1]
    probTrans[i, 1] ~ dunif(0, 1)
    probTrans[i, 2] ~ dunif(0, 1 - probTrans[i, 1])
    probTrans[i, 3] <- 1 - probTrans[i, 1] - probTrans[i, 2]
  }

})
Rmodel <- nimbleModel(nc, data = list(x = x),
                 inits = list(
                   init = init,
                   probObs = probObs,
                   probTrans = probTrans
                 ),
                 buildDerivs=TRUE)
Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)
Cmodel$calculate()

nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c(Rmodel$expandNodeNames('init[1:3]'),
                                                                       Rmodel$expandNodeNames('probObs[1:3, 1:2]'),
                                                                       Rmodel$expandNodeNames('probTrans[1:3, 1:3]')
                                                                       ))
v1_case1 <- list(arg1 = c(init[1:3],  probObs[1:3, 1:2], probTrans[1:3, 1:3]))
v2_case1 <- list(arg1 = c(init2[1:3], probObs2[1:3, 1:2],  probTrans2[1:3, 1:3]))

model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                          nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                          order = 0:2)
#######

######################
#### dHMMo case ####

x <- c(1, 1, 1, 2, 2)

init <- c(0.4, 0.2, 0.4)
probObs <- array(
       c(0.95, 0.05, 0.8, 0.05, 0.95, 0.2,
         0.95, 0.05, 0.8, 0.05, 0.95, 0.2,
         0.9, 0.05, 0.8, 0.1, 0.95, 0.2,
         0.95, 0.1, 0.8, 0.05, 0.9, 0.2,
         0.95, 0.05, 0.7, 0.05, 0.95, 0.3),
       c(3, 2, 5))

probTrans <- t(array(
          c(0.3, 0.4, 0.2,
            0.1, 0.1, 0.8,
            0.05, 0.05, 0.9),
          c(3,3)))

init2 <- c(0.6, 0.1, 0.3)
probObs2 <- array(
       c(0.6, 0.05, 0.8, 0.4, 0.95, 0.2,
         0.8, 0.05, 0.6, 0.2, 0.95, 0.4,
         0.9, 0.05, 0.8, 0.1, 0.95, 0.2,
         0.9, 0.1, 0.8, 0.1, 0.9, 0.2,
         0.95, 0.05, 0.4, 0.05, 0.95, 0.6),
       c(3, 2, 5))

probTrans2 <- t(array(
          c(0.4, 0.4, 0.2,
            0.05, 0.25, 0.7,
            0.05, 0.15, 0.8),
          c(3,3)))

nc <- nimbleCode({
  x[1:5] ~ dHMMo(init[1:3], probObs = probObs[1:3,1:2,1:5],
                  probTrans = probTrans[1:3, 1:3], len = 5, checkRowSums = 0)
  for (i in 1:2) {
    init[i] ~ dunif(0, 1)
  }
  init[3] <- 1 - init[1] - init[2]

  for (i in 1:3) {
    for (j in 1:5) {
      probObs[i, 1, j] ~ dunif(0, 1)
      probObs[i, 2, j] <- 1 - probObs[i, 1, j]
    }
    probTrans[i, 1] ~ dunif(0, 1)
    probTrans[i, 2] ~ dunif(0, 1 - probTrans[i, 1])
    probTrans[i, 3] <- 1 - probTrans[i, 1] - probTrans[i, 2]
  }

})
Rmodel <- nimbleModel(nc, data = list(x = x),
                 inits = list(
                   init = init,
                   probObs = probObs,
                   probTrans = probTrans
                 ),
                 buildDerivs=TRUE)
Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)
Cmodel$calculate()

nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c(Rmodel$expandNodeNames('init[1:3]'),
                                                                       Rmodel$expandNodeNames('probObs[1:3, 1:2, 1:5]'),
                                                                       Rmodel$expandNodeNames('probTrans[1:3, 1:3]')
                                                                       ))
v1_case1 <- list(arg1 = c(init[1:3],  probObs[1:3, 1:2, 1:5], probTrans[1:3, 1:3]))
v2_case1 <- list(arg1 = c(init2[1:3], probObs2[1:3, 1:2, 1:5],  probTrans2[1:3, 1:3]))

model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                          nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                          order = 0:2)
#######


######################
#### dDHMM case ####

x <- c(1, 1, 1, 2, 2)

init <- c(0.4, 0.2, 0.4)
probObs <- t(array(
         c(0.9, 0.1,
           0.1, 0.9,
           0.8, 0.2),
         c(2, 3)))

probTrans <- array(
        c(0.6, 0.05, 0.05, 0.3, 0.65, 0.25, 0.1, 0.3, 0.7,
          0.6, 0.05, 0.05, 0.2, 0.65, 0.05, 0.2, 0.3, 0.9,
          0.6, 0.05, 0.05, 0.2, 0.65, 0.05, 0.2, 0.3, 0.9,
          0.6, 0.05, 0.20, 0.3, 0.65, 0.05, 0.1, 0.3, 0.75
          ),
        c(3,3,4))

init2 <- c(0.6, 0.1, 0.3)
probObs2 <- t(array(
         c(0.9, 0.1,
           0.1, 0.9,
           0.7, 0.3),
         c(2, 3)))
probTrans2 <- array(
        c(0.5, 0.05, 0.02, 0.4, 0.65, 0.28, 0.1, 0.3, 0.7,
          0.5, 0.05, 0.02, 0.3, 0.75, 0.08, 0.2, 0.2, 0.9,
          0.6, 0.05, 0.05, 0.2, 0.75, 0.05, 0.2, 0.2, 0.9,
          0.6, 0.05, 0.20, 0.3, 0.65, 0.05, 0.1, 0.3, 0.75
          ),
        c(3,3,4))

nc <- nimbleCode({
  x[1:5] ~ dDHMM(init[1:3], probObs = probObs[1:3,1:2],
                  probTrans = probTrans[1:3, 1:3, 1:4], len = 5, checkRowSums = 0)
  for (i in 1:2) {
    init[i] ~ dunif(0, 1)
  }
  init[3] <- 1 - init[1] - init[2]

  for (i in 1:3) {
    probObs[i, 1] ~ dunif(0, 1)
    probObs[i, 2] <- 1 - probObs[i, 1]
    for (k in 1:4) {
      probTrans[i, 1, k] ~ dunif(0, 1)
      probTrans[i, 2, k] ~ dunif(0, 1 - probTrans[i, 1, k])
      probTrans[i, 3, k] <- 1 - probTrans[i, 1, k] - probTrans[i, 2, k]
    }
  }

})
Rmodel <- nimbleModel(nc, data = list(x = x),
                 inits = list(
                   init = init,
                   probObs = probObs,
                   probTrans = probTrans
                 ),
                 buildDerivs=TRUE)
Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)
Cmodel$calculate()

nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c(Rmodel$expandNodeNames('init[1:3]'),
                                                                       Rmodel$expandNodeNames('probObs[1:3, 1:2, 1:5]'),
                                                                       Rmodel$expandNodeNames('probTrans[1:3, 1:3, 1:4]')
                                                                       ))
v1_case1 <- list(arg1 = c(init[1:3],  probObs[1:3, 1:2, 1:5], probTrans[1:3, 1:3, 1:4]))
v2_case1 <- list(arg1 = c(init2[1:3], probObs2[1:3, 1:2, 1:5],  probTrans2[1:3, 1:3, 1:4]))

model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                          nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                          order = 0:2)
#######


######################
#### dDHMMo case ####

x <- c(1, 1, 1, 2, 2)

init <- c(0.4, 0.2, 0.4)
probObs <- array(
       c(0.95, 0.05, 0.8, 0.05, 0.95, 0.2,
         0.95, 0.05, 0.8, 0.05, 0.95, 0.2,
         0.9, 0.05, 0.8, 0.1, 0.95, 0.2,
         0.95, 0.1, 0.8, 0.05, 0.9, 0.2,
         0.95, 0.05, 0.7, 0.05, 0.95, 0.3),
       c(3, 2, 5))

probTrans <- array(
        c(0.6, 0.05, 0.05, 0.3, 0.65, 0.25, 0.1, 0.3, 0.7,
          0.6, 0.05, 0.05, 0.2, 0.65, 0.05, 0.2, 0.3, 0.9,
          0.6, 0.05, 0.05, 0.2, 0.65, 0.05, 0.2, 0.3, 0.9,
          0.6, 0.05, 0.20, 0.3, 0.65, 0.05, 0.1, 0.3, 0.75
          ),
        c(3,3,4))

init2 <- c(0.6, 0.1, 0.3)
probObs2 <- array(
       c(0.6, 0.05, 0.8, 0.4, 0.95, 0.2,
         0.8, 0.05, 0.6, 0.2, 0.95, 0.4,
         0.9, 0.05, 0.8, 0.1, 0.95, 0.2,
         0.9, 0.1, 0.8, 0.1, 0.9, 0.2,
         0.95, 0.05, 0.4, 0.05, 0.95, 0.6),
       c(3, 2, 5))

probTrans2 <- array(
        c(0.5, 0.05, 0.02, 0.4, 0.65, 0.28, 0.1, 0.3, 0.7,
          0.5, 0.05, 0.02, 0.3, 0.75, 0.08, 0.2, 0.2, 0.9,
          0.6, 0.05, 0.05, 0.2, 0.75, 0.05, 0.2, 0.2, 0.9,
          0.6, 0.05, 0.20, 0.3, 0.65, 0.05, 0.1, 0.3, 0.75
          ),
        c(3,3,4))

nc <- nimbleCode({
  x[1:5] ~ dDHMMo(init[1:3], probObs = probObs[1:3,1:2,1:5],
                  probTrans = probTrans[1:3, 1:3, 1:4], len = 5, checkRowSums = 0)
  for (i in 1:2) {
    init[i] ~ dunif(0, 1)
  }
  init[3] <- 1 - init[1] - init[2]

  for (i in 1:3) {
    for (j in 1:5) {
      probObs[i, 1, j] ~ dunif(0, 1)
      probObs[i, 2, j] <- 1 - probObs[i, 1, j]
    }
    for (k in 1:4) {
      probTrans[i, 1, k] ~ dunif(0, 1)
      probTrans[i, 2, k] ~ dunif(0, 1 - probTrans[i, 1, k])
      probTrans[i, 3, k] <- 1 - probTrans[i, 1, k] - probTrans[i, 2, k]
    }
  }

})
Rmodel <- nimbleModel(nc, data = list(x = x),
                 inits = list(
                   init = init,
                   probObs = probObs,
                   probTrans = probTrans
                 ),
                 buildDerivs=TRUE)
Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)
Cmodel$calculate()

nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c(Rmodel$expandNodeNames('init[1:3]'),
                                                                       Rmodel$expandNodeNames('probObs[1:3, 1:2, 1:5]'),
                                                                       Rmodel$expandNodeNames('probTrans[1:3, 1:3, 1:4]')
                                                                       ))
v1_case1 <- list(arg1 = c(init[1:3],  probObs[1:3, 1:2, 1:5], probTrans[1:3, 1:3, 1:4]))
v2_case1 <- list(arg1 = c(init2[1:3], probObs2[1:3, 1:2, 1:5],  probTrans2[1:3, 1:3, 1:4]))

model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                          nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                          order = 0:2)
#######

######################
#### dDynOcc_vvm case ####

x <- matrix(c(0,0,NA,0,
              1,1,1,0,
              0,0,0,0,
              0,0,1,0,
              0,0,0,NA), nrow = 4)
start <- c(1,1,2,1)
end <- c(5,5,5,4)

init <- 0.7
probPersist <- c(0.4, 0.4, 0.1)
probColonize <- c(0.4, 0.2, 0.1)
p <- matrix(rep(c(0.8, 0.7, 0.8, 0.8, 0.9), each = 4), nrow = 4, byrow =TRUE)

init2 <- 0.9
probPersist2 <- c(0.4, 0.4, 0.1)
probColonize2 <- c(0.4, 0.2, 0.1)
p2 <- matrix(rep(c(0.7, 0.5, 0.3, 0.8, 0.66), each = 4), nrow = 4, byrow =TRUE)


nc <- nimbleCode({
    x[1:4, 1:5] ~ dDynOcc_vvm(init,
                              probPersist[1:3],
                              probColonize[1:3],
                              p[1:4,1:5],
                              start[1:4], end[1:4])

  init ~ dunif(0, 1)
  for (i in 1:3) {
    probColonize[i] ~ dunif(0, 1)
    probPersist[i] ~ dunif(0, 1)
  }
  for (i in 1:4) {
    for (j in 1:5) {
      p[i, j] ~ dunif(0, 1)
    }
  }

})

Rmodel <- nimbleModel(nc, data = list(x = x),
                      constants = list(start = start, end = end),
                 inits = list(
                   init = init,
                   p = p,
                   probColonize = probColonize,
                   probPersist = probPersist
                 ),
                 buildDerivs=TRUE)
Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)
Cmodel$calculate()

nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c(
       "init",
       Rmodel$expandNodeNames('p[1:4, 1:5]'),
       Rmodel$expandNodeNames('probColonize[1:3]'),
       Rmodel$expandNodeNames('probPersist[1:3]')
   ))
v1_case1 <- list(arg1 = c(init, p, probColonize, probPersist))
v2_case1 <- list(arg1 = c(init, p2, probColonize2, probPersist2))

model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                          nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                          order = 0:2)



######################
#### dDynOcc_vvv case ####

x <- matrix(c(0,0,NA,0,
              1,1,1,0,
              0,0,0,0,
              0,0,1,0,
              0,0,0,NA), nrow = 4)
start <- c(1,1,2,1)
end <- c(5,5,5,4)

init <- 0.7
probPersist <- c(0.4, 0.4, 0.1)
probColonize <- c(0.4, 0.2, 0.1)
p <- c(0.8, 0.7, 0.8, 0.8)

init2 <- 0.9
probPersist2 <- c(0.4, 0.4, 0.1)
probColonize2 <- c(0.4, 0.2, 0.1)
p2 <- c(0.7, 0.5, 0.3, 0.8)


nc <- nimbleCode({
    x[1:4, 1:5] ~ dDynOcc_vvv(init,
                              probPersist[1:3],
                              probColonize[1:3],
                              p[1:4],
                              start[1:4], end[1:4])

  init ~ dunif(0, 1)
  for (i in 1:3) {
    probColonize[i] ~ dunif(0, 1)
    probPersist[i] ~ dunif(0, 1)
  }
  for (i in 1:4) {
    p[i] ~ dunif(0, 1)
  }

})

Rmodel <- nimbleModel(nc, data = list(x = x),
                      constants = list(start = start, end = end),
                 inits = list(
                   init = init,
                   p = p,
                   probColonize = probColonize,
                   probPersist = probPersist
                 ),
                 buildDerivs=TRUE)
Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)
Cmodel$calculate()

nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c(
       "init",
       Rmodel$expandNodeNames('p[1:4]'),
       Rmodel$expandNodeNames('probColonize[1:3]'),
       Rmodel$expandNodeNames('probPersist[1:3]')
   ))
v1_case1 <- list(arg1 = c(init, p, probColonize, probPersist))
v2_case1 <- list(arg1 = c(init, p2, probColonize2, probPersist2))

model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                          nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                          order = 0:2)

######################
#### dDynOcc_vvs case ####

x <- matrix(c(0,0,NA,0,
              1,1,1,0,
              0,0,0,0,
              0,0,1,0,
              0,0,0,NA), nrow = 4)
start <- c(1,1,2,1)
end <- c(5,5,5,4)

init <- 0.7
probPersist <- c(0.4, 0.4, 0.1)
probColonize <- c(0.4, 0.2, 0.1)
p <- c(0.8)

init2 <- 0.9
probPersist2 <- c(0.4, 0.4, 0.1)
probColonize2 <- c(0.4, 0.2, 0.1)
p2 <- c(0.7)


nc <- nimbleCode({
    x[1:4, 1:5] ~ dDynOcc_vvs(init,
                              probPersist[1:3],
                              probColonize[1:3],
                              p,
                              start[1:4], end[1:4])

  init ~ dunif(0, 1)
  for (i in 1:3) {
    probColonize[i] ~ dunif(0, 1)
    probPersist[i] ~ dunif(0, 1)
  }
  p ~ dunif(0, 1)

})

Rmodel <- nimbleModel(nc, data = list(x = x),
                      constants = list(start = start, end = end),
                 inits = list(
                   init = init,
                   p = p,
                   probColonize = probColonize,
                   probPersist = probPersist
                 ),
                 buildDerivs=TRUE)
Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)
Cmodel$calculate()

nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c(
       "init", "p",
       Rmodel$expandNodeNames('probColonize[1:3]'),
       Rmodel$expandNodeNames('probPersist[1:3]')
   ))
v1_case1 <- list(arg1 = c(init, p, probColonize, probPersist))
v2_case1 <- list(arg1 = c(init, p2, probColonize2, probPersist2))

model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                          nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                          order = 0:2)



# reset options before finishing
nimbleOptions(enableDerivs = EDopt)
nimbleOptions(buildModelDerivs = BMDopt)





#### Notes:
#' dNmixture appears not to work. When the model is first defined I get the following message:
#'     [Note] Detected use of function(s) that are not supported for derivative tracking in a function or method for which buildDerivs has been requested: qpois.
#'     Then model compilation fails.
#' in dCJS_*v variations, I had to manually specify to only do derivs for probCapture[2:n] since element probCapture[1] is
#'     ignored in the likelihood. Works fine with this change
#' hit issues with d*HMM*. The derivatives work fine, except if probabilities equal to 0 are included in the
#'     transition matrices. Even then the error only occurs during model_calculate_test_case.
#'     The error I get is: "Error in if (!all_result) { : missing value where TRUE/FALSE needed"
#'     which I beleive stems from an NA derivative. (Could be helpful to have a more informative error msg here.)
#'
#'Got the following error in dDHMMo test, and an equivalent error in dDynOcc_vvm, dDynOcc_vvv
#'Detected some values out of relative tolerance  (RC order 0) :  as.numeric(first)   as.numeric(others[[i]]) .
# [1] 2.140095e-01 2.140095e-01 2.075091e-15
# ******************
# Some C-to-R derivatives to not match for order 0Called from: test_AD2_oneCall(Rfxn, Cfxn, recordArgs = v1, testArgs = v2,
#     order = order, Rmodel = Rmodel, Cmodel = Cmodel, recordInits = varValues,
#     testInits = varValues2, nodesToChange = c(nodesList$updateNodes),
#     ...)
#'
#'
