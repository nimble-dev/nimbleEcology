# Testing examples:

# install nimble from CRAN:
#devtools::install_cran("nimble", force = TRUE)
# devtools::install_github("nimble-dev/nimble", ref = "devel", subdir = "packages/nimble")
# install nimbleEcology from branch Nmixture-AD: devtools::install_github("nimble-dev/nimbleEcology", ref = "Nmixture-AD")

# load nimble's testing tools
#library(nimble)
library(nimbleEcology)
source(system.file('test_utils.R', package = 'nimbleEcology'))
source(system.file('AD_test_utils.R', package = 'nimbleEcology'))
# source("../nimble/packages/nimble/tests/testthat/test_utils.R")
# source("../nimble/packages/nimble/tests/testthat/AD_test_utils.R")

EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)

test_that("dOcc works with AD",
{
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

  nodesList_case1 <-
    setup_update_and_constant_nodes_for_tests(Rmodel, c('probOcc', 'probDetect'))
  v1_case1 <- list(arg1 = c(0.6, 0.4)) # taping values for probOcc and probDetect
  v2_case1 <- list(arg1 = c(0.65, 0.35)) # testing values for probOcc and probDetect
  RCrelTol = c(1e-15, 1e-8, 1e-3, 1e-14)

  res <- model_calculate_test_case(Rmodel, Cmodel,
                            model_calculate_test, nodesList_case1,
                            v1_case1, v2_case1,
                            0:2) # lots of output numbers with no warning messages means it passes.

  # Missing values
  dat2 <- c(1,NA,0,0) # A vector of observations

  nc <- nimbleCode({
    x[1:4] ~ dOcc_s(probOcc, probDetect, len = 4)
    probOcc ~ dunif(0,1)
    probDetect ~ dunif(0,1)
  })

  Rmodel_na <- nimbleModel(nc, data = list(x = dat2),
                        inits = list(probOcc = probOcc,
                                     probDetect = probDetect),
                        buildDerivs=TRUE)

  Cmodel_na <- compileNimble(Rmodel_na)

  nodesList_case_na <-
    setup_update_and_constant_nodes_for_tests(Rmodel_na, c('probOcc', 'probDetect'))
  v1_case1 <- list(arg1 = c(0.6, 0.4)) # taping values for probOcc and probDetect
  v2_case1 <- list(arg1 = c(0.65, 0.35)) # testing values for probOcc and probDetect
  RCrelTol = c(1e-15, 1e-8, 1e-3, 1e-14)

  res_na <- model_calculate_test_case(Rmodel_na, Cmodel_na,
                            model_calculate_test, nodesList_case_na,
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

  res <- model_calculate_test_case(Rmodel, Cmodel,
                            model_calculate_test, nodesList_case1,
                            v1_case1, v2_case1,
                            0:2)

  # Missing values
  x2 <- c(1,0,NA,1,0)
  Rmodel_na <- nimbleModel(nc, data = list(x = x2),
                        inits = list(probOcc = probOcc,
                                     probDetect = probDetect),
                        buildDerivs=TRUE)
  Rmodel_na$calculate()

  Cmodel_na <- compileNimble(Rmodel_na)

  nodesList_case_na <- setup_update_and_constant_nodes_for_tests(Rmodel_na,
                                                               c('probOcc', Rmodel_na$expandNodeNames('probDetect[1:5]')))
  v1_case1 <- list(arg1 = c(probOcc, probDetect)) # taping values for probOcc and probDetect
  v2_case1 <- list(arg1 = c(probOcc2, probDetect2)) # testing values for probOcc and probDetect

  res <- model_calculate_test_case(Rmodel_na, Cmodel_na,
                            model_calculate_test, nodesList_case_na,
                            v1_case1, v2_case1,
                            0:2)
})

test_that ("dNmixture works with AD", {
##########################
#### dNmixture_s case ####

  x <- c(7, 7, 6, 9, 10)
  lambda <- 15
  prob <- 0.7

  lambda2 <- 18
  prob2 <- 0.5


  nc <- nimbleCode({
    x[1:5] ~ dNmixtureAD_s(lambda, prob,
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
  res <- model_calculate_test_case(Rmodel, Cmodel,
                                   model_calculate_test, nodesList_case1,
                                   v1_case1, v2_case1,
                                   0:2, RCrelTol = c(2e-15, 1e-8, 1e-3, 1e-14))

  # Missing values
  xna <- c(7, 7, NA, 9, 10)
  Rmodel_na <- nimbleModel(nc, data = list(x = xna),
                        inits = list(prob = prob,
                                     lambda = lambda),
                        buildDerivs=TRUE)
  Rmodel_na$calculate()

  Cmodel_na <- compileNimble(Rmodel_na)
  Cmodel_na$calculate()

  nodesList_case1_na <- setup_update_and_constant_nodes_for_tests(Rmodel_na, c('prob', 'lambda'))

  res_na <- model_calculate_test_case(Rmodel_na, Cmodel_na,
                            model_calculate_test, nodesList_case1_na,
                            v1_case1, v2_case1,
                            0:2, RCrelTol = c(2e-15, 1e-8, 1e-3, 1e-14))

##########################
#### dNmixture_BNB_s case ####

  x <- c(7, 7, 6, 9, 10)
  lambda <- 15
  prob <- 0.7
  theta <- 1.1

  lambda2 <- 18
  prob2 <- 0.5
  theta2 <- 1.2

  nc <- nimbleCode({
    x[1:5] ~ dNmixtureAD_BNB_s(lambda, prob, theta = theta,
                             Nmin = 0, Nmax = 100, len = 5)
    prob ~ dunif(0, 1)
    lambda ~ dunif(0, 100)
  })

  Rmodel <- nimbleModel(nc, data = list(x = x),
                        inits = list(prob = prob,
                                     lambda = lambda,
                                     theta = theta),
                        buildDerivs=TRUE)
  Rmodel$calculate()

  Cmodel <- compileNimble(Rmodel)
  Cmodel$calculate()

  nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c('prob', 'lambda', 'theta'))
  v1_case1 <- list(arg1 = c(prob, lambda, theta)) # taping values for prob and lambda
  v2_case1 <- list(arg1 = c(prob2, lambda2, theta2)) # testing values for prob and lambda

  res <- model_calculate_test_case(Rmodel, Cmodel,
                            model_calculate_test, nodesList_case1,
                            v1_case1, v2_case1,
                            0:2, RCrelTol = c(2e-15, 1e-8, 1e-3, 1e-14))

  # NA handling
  x <- c(7, 7, NA, 9, 10)
  nc <- nimbleCode({
    x[1:5] ~ dNmixtureAD_BNB_s(lambda, prob, theta = theta,
                             Nmin = 0, Nmax = 100, len = 5)
    prob ~ dunif(0, 1)
    lambda ~ dunif(0, 100)
  })

  Rmodel <- nimbleModel(nc, data = list(x = x),
                        inits = list(prob = prob,
                                     lambda = lambda,
                                     theta = theta),
                        buildDerivs=TRUE)
  Rmodel$calculate()

  Cmodel <- compileNimble(Rmodel)
  Cmodel$calculate()

  nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c('prob', 'lambda', 'theta'))
  v1_case1 <- list(arg1 = c(prob, lambda, theta)) # taping values for prob and lambda
  v2_case1 <- list(arg1 = c(prob2, lambda2, theta2)) # testing values for prob and lambda

  res <- model_calculate_test_case(Rmodel, Cmodel,
                            model_calculate_test, nodesList_case1,
                            v1_case1, v2_case1,
                            0:2, RCrelTol = c(2e-15, 1e-8, 1e-3, 1e-14))

##############################
#### dNmixture_BBP_s case ####

  x <- c(7, 7, 6, 9, 10)
  lambda <- 15
  prob <- 0.7
  s <- 0.93

  lambda2 <- 18
  prob2 <- 0.5
  s2 <- 1.111

  nc <- nimbleCode({
    x[1:5] ~ dNmixtureAD_BBP_s(lambda, prob, s = s,
                               Nmin = 0, Nmax = 100, len = 5)
    prob ~ dunif(0, 1)
    lambda ~ dunif(0, 100)
  })

  Rmodel <- nimbleModel(nc, data = list(x = x),
                        inits = list(prob = prob,
                                     lambda = lambda,
                                     s = s),
                        buildDerivs = TRUE)
  Rmodel$calculate()

  Cmodel <- compileNimble(Rmodel)
  Cmodel$calculate()

  nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c('prob', 'lambda', 's'))
  v1_case1 <- list(arg1 = c(prob, lambda, s)) # taping values for prob and lambda
  v2_case1 <- list(arg1 = c(prob2, lambda2, s2)) # testing values for prob and lambda

  res <- model_calculate_test_case(Rmodel, Cmodel,
                            model_calculate_test, nodesList_case1,
                            v1_case1, v2_case1,
                            0:2, RCrelTol = c(2e-15, 1e-8, 1e-3, 1e-14))

  # Missing value handling
  x <- c(7, 7, NA, 9, 10)
  Rmodel <- nimbleModel(nc, data = list(x = x),
                        inits = list(prob = prob,
                                     lambda = lambda,
                                     s = s),
                        buildDerivs = TRUE)
  Rmodel$calculate()

  Cmodel <- compileNimble(Rmodel)
  Cmodel$calculate()

  nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c('prob', 'lambda', 's'))
  v1_case1 <- list(arg1 = c(prob, lambda, s)) # taping values for prob and lambda
  v2_case1 <- list(arg1 = c(prob2, lambda2, s2)) # testing values for prob and lambda

  res <- model_calculate_test_case(Rmodel, Cmodel,
                            model_calculate_test, nodesList_case1,
                            v1_case1, v2_case1,
                            0:2, RCrelTol = c(2e-15, 1e-8, 1e-3, 1e-14))

##############################
#### dNmixture_BBNB_s case ####

  x <- c(7, 7, 6, 9, 10)
  lambda <- 15
  prob <- 0.7
  theta <- 1.1
  s <- 0.93

  lambda2 <- 18
  prob2 <- 0.5
  theta2 <- 1.2
  s2 <- 1.111

  nc <- nimbleCode({
    x[1:5] ~ dNmixtureAD_BBNB_s(lambda, prob, theta = theta, s = s,
                              Nmin = 0, Nmax = 100, len = 5)
    prob ~ dunif(0, 1)
    lambda ~ dunif(0, 100)
  })

  Rmodel <- nimbleModel(nc, data = list(x = x),
                        inits = list(prob = prob,
                                     lambda = lambda,
                                     theta = theta, s = s),
                        buildDerivs=TRUE)
  Rmodel$calculate()

  Cmodel <- compileNimble(Rmodel)
  Cmodel$calculate()

  nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c('prob', 'lambda', 'theta', 's'))
  v1_case1 <- list(arg1 = c(prob, lambda, theta, s)) # taping values for prob and lambda
  v2_case1 <- list(arg1 = c(prob2, lambda2, theta2, s2)) # testing values for prob and lambda

  res <- model_calculate_test_case(Rmodel, Cmodel,
                            model_calculate_test, nodesList_case1,
                            v1_case1, v2_case1,
                            0:2, RCrelTol = c(2e-15, 1e-8, 1e-3, 1e-14))

##########################
#### dNmixture_v case ####

  x <- c(7, 7, 6, 9, 10)
  lambda <- 15
  prob <- c(0.6, 0.6, 0.4, 0.9, 0.8)

  lambda2 <- 18
  prob2 <-  c(0.65, 0.65, 0.45, 0.95, 0.85)


  nc <- nimbleCode({
    x[1:5] ~ dNmixtureAD_v(lambda, prob[1:5],
                           Nmin = 0, Nmax = 100, len = 5)
    for (i in 1:5) {
      prob[i] ~ dunif(0, 1)
    }
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

  res <- model_calculate_test_case(Rmodel, Cmodel,
                            model_calculate_test, nodesList_case1,
                            v1_case1, v2_case1,
                            0:2, RCrelTol = c(2e-15, 1e-8, 1e-3, 1e-14))

  # Missing values
  xna <- c(7, 7, NA, 9, 10)
  Rmodel_na <- nimbleModel(nc, data = list(x = xna),
                        inits = list(prob = prob,
                                     lambda = lambda),
                        buildDerivs=TRUE)
  Rmodel_na$calculate()

  Cmodel_na <- compileNimble(Rmodel_na)
  Cmodel_na$calculate()

  nodesList_case1_na <- setup_update_and_constant_nodes_for_tests(Rmodel_na, c('prob', 'lambda'))

  res_na <- model_calculate_test_case(Rmodel_na, Cmodel_na,
                            model_calculate_test, nodesList_case1_na,
                            v1_case1, v2_case1,
                            0:2, RCrelTol = c(2e-15, 1e-8, 1e-3, 1e-14))

##############################
#### dNmixture_BNB_v case ####

  x <- c(7, 7, 6, 9, 10)
  lambda <- 15
  prob <- c(0.6, 0.6, 0.4, 0.9, 0.8)
  theta = 1.1


  lambda2 <- 18
  prob2 <-  c(0.65, 0.65, 0.45, 0.95, 0.85)
  theta2 <- 1.2


  nc <- nimbleCode({
    x[1:5] ~ dNmixtureAD_BNB_v(lambda, prob[1:5], theta = theta,
                               Nmin = 0, Nmax = 100, len = 5)
    for (i in 1:5) {
      prob[i] ~ dunif(0, 1)
    }
    lambda ~ dunif(0, 100)
  })

  Rmodel <- nimbleModel(nc, data = list(x = x),
                        inits = list(prob = prob,
                                     lambda = lambda,
                                     theta = theta),
                        buildDerivs=TRUE)
  Rmodel$calculate()

  Cmodel <- compileNimble(Rmodel)
  Cmodel$calculate()

  nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c('prob', 'lambda', 'theta'))
  v1_case1 <- list(arg1 = c(prob, lambda, theta)) # taping values for prob and lambda
  v2_case1 <- list(arg1 = c(prob2, lambda2, theta2)) # testing values for prob and lambda

  res <- model_calculate_test_case(Rmodel, Cmodel,
                            model_calculate_test, nodesList_case1,
                            v1_case1, v2_case1,
                            0:2, RCrelTol = c(2e-15, 1e-8, 1e-3, 1e-14))

  # NA handling
  x <- c(7, 7, NA, 9, 10)
  Rmodel <- nimbleModel(nc, data = list(x = x),
                        inits = list(prob = prob,
                                     lambda = lambda,
                                     theta = theta),
                        buildDerivs=TRUE)
  Rmodel$calculate()

  Cmodel <- compileNimble(Rmodel)
  Cmodel$calculate()

  nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c('prob', 'lambda', 'theta'))
  v1_case1 <- list(arg1 = c(prob, lambda, theta)) # taping values for prob and lambda
  v2_case1 <- list(arg1 = c(prob2, lambda2, theta2)) # testing values for prob and lambda

  res <- model_calculate_test_case(Rmodel, Cmodel,
                            model_calculate_test, nodesList_case1,
                            v1_case1, v2_case1,
                            0:2, RCrelTol = c(2e-15, 1e-8, 1e-3, 1e-14))

##############################
#### dNmixture_BBP_v case ####

  x <- c(7, 7, 6, 9, 10)
  lambda <- 15
  prob <- c(0.6, 0.6, 0.4, 0.9, 0.8)
  s <- 0.9


  lambda2 <- 18
  prob2 <-  c(0.65, 0.65, 0.45, 0.95, 0.85)
  s2 <- 1.2


  nc <- nimbleCode({
    x[1:5] ~ dNmixtureAD_BBP_v(lambda, prob[1:5], s = s,
                               Nmin = 0, Nmax = 100, len = 5)
    for (i in 1:5) {
      prob[i] ~ dunif(0, 1)
    }
    lambda ~ dunif(0, 100)
  })

  Rmodel <- nimbleModel(nc, data = list(x = x),
                        inits = list(prob = prob,
                                     lambda = lambda,
                                     s = s),
                        buildDerivs=TRUE)
  Rmodel$calculate()

  Cmodel <- compileNimble(Rmodel)
  Cmodel$calculate()

  nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c('prob', 'lambda', 's'))
  v1_case1 <- list(arg1 = c(prob, lambda, s)) # taping values for prob and lambda
  v2_case1 <- list(arg1 = c(prob2, lambda2, s2)) # testing values for prob and lambda

  res <- model_calculate_test_case(Rmodel, Cmodel,
                            model_calculate_test, nodesList_case1,
                            v1_case1, v2_case1,
                            0:2, RCrelTol = c(2e-15, 1e-8, 1e-3, 1e-14))

  # Missing value handling
  x <- c(7, 7, NA, 9, 10)
  Rmodel <- nimbleModel(nc, data = list(x = x),
                        inits = list(prob = prob,
                                     lambda = lambda,
                                     s = s),
                        buildDerivs=TRUE)
  Rmodel$calculate()

  Cmodel <- compileNimble(Rmodel)
  Cmodel$calculate()

  nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c('prob', 'lambda', 's'))
  v1_case1 <- list(arg1 = c(prob, lambda, s)) # taping values for prob and lambda
  v2_case1 <- list(arg1 = c(prob2, lambda2, s2)) # testing values for prob and lambda

  res <- model_calculate_test_case(Rmodel, Cmodel,
                            model_calculate_test, nodesList_case1,
                            v1_case1, v2_case1,
                            0:2, RCrelTol = c(2e-15, 1e-8, 1e-3, 1e-14))

##############################
#### dNmixture_BBNB_v case ####

  x <- c(7, 7, 6, 9, 10)
  lambda <- 15
  prob <- c(0.6, 0.6, 0.4, 0.9, 0.8)
  theta <- 1.1
  s <- 0.9


  lambda2 <- 18
  prob2 <-  c(0.65, 0.65, 0.45, 0.95, 0.85)
  theta2 <- 1.2
  s2 <- 1.2


  nc <- nimbleCode({
    x[1:5] ~ dNmixtureAD_BBNB_v(lambda, prob[1:5], theta = theta, s = s,
                              Nmin = 0, Nmax = 100, len = 5)
    for (i in 1:5) {
      prob[i] ~ dunif(0, 1)
    }
    lambda ~ dunif(0, 100)
  })

  Rmodel <- nimbleModel(nc, data = list(x = x),
                        inits = list(prob = prob,
                                     lambda = lambda,
                                     theta = theta, s = s),
                        buildDerivs=TRUE)
  Rmodel$calculate()

  Cmodel <- compileNimble(Rmodel)
  Cmodel$calculate()

  nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c('prob', 'lambda', 'theta', 's'))
  v1_case1 <- list(arg1 = c(prob, lambda, theta, s)) # taping values for prob and lambda
  v2_case1 <- list(arg1 = c(prob2, lambda2, theta2, s2)) # testing values for prob and lambda

  res <- model_calculate_test_case(Rmodel, Cmodel,
                            model_calculate_test, nodesList_case1,
                            v1_case1, v2_case1,
                            0:2, RCrelTol = c(2e-15, 1e-8, 1e-3, 1e-14))

##########################
#### dNmixture_BNB_oneObs case ####

  x <- 8
  lambda <- 15
  prob <- 0.7
  theta <- 1.1

  lambda2 <- 18
  prob2 <- 0.66
  theta2 <- 1.4

  nc <- nimbleCode({
    x ~ dNmixtureAD_BNB_oneObs(lambda, prob, theta = theta,
                             Nmin = 0, Nmax = 100)
    prob ~ dunif(0, 1)
    lambda ~ dunif(0, 100)
  })

  Rmodel <- nimbleModel(nc, data = list(x = x),
                        inits = list(prob = prob, theta = theta,
                                     lambda = lambda),
                        buildDerivs=TRUE)
  Rmodel$calculate()

  Cmodel <- compileNimble(Rmodel)
  Cmodel$calculate()

  nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c('prob', 'lambda', 'theta'))
  v1_case1 <- list(arg1 = c(prob, lambda, theta)) # taping values for prob and lambda
  v2_case1 <- list(arg1 = c(prob2, lambda2, theta2)) # testing values for prob and lambda

  res <- model_calculate_test_case(Rmodel, Cmodel,
                            model_calculate_test, nodesList_case1,
                            v1_case1, v2_case1,
                            0:2, RCrelTol = c(2e-15, 1e-8, 1e-3, 1e-14))

##########################
#### dNmixture_BBP_oneObs case ####

  x <- 8
  lambda <- 15
  prob <- 0.7
  s <- 0.8

  lambda2 <- 18
  prob2 <- 0.66
  s <- 1.1

  nc <- nimbleCode({
    x ~ dNmixtureAD_BBP_oneObs(lambda, prob, s = s,
                             Nmin = 0, Nmax = 100)
    prob ~ dunif(0, 1)
    lambda ~ dunif(0, 100)
  })

  Rmodel <- nimbleModel(nc, data = list(x = x),
                        inits = list(prob = prob, s=s,
                                     lambda = lambda),
                        buildDerivs=TRUE)
  Rmodel$calculate()

  Cmodel <- compileNimble(Rmodel)
  Cmodel$calculate()

  nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c('prob', 'lambda', 's'))
  v1_case1 <- list(arg1 = c(prob, lambda, s2)) # taping values for prob and lambda
  v2_case1 <- list(arg1 = c(prob2, lambda2, s2)) # testing values for prob and lambda

  res <- model_calculate_test_case(Rmodel, Cmodel,
                            model_calculate_test, nodesList_case1,
                            v1_case1, v2_case1,
                            0:2, RCrelTol = c(2e-15, 1e-8, 1e-3, 1e-14))

##########################
#### dNmixture_BBNB_oneObs case ####

  x <- 8
  lambda <- 15
  prob <- 0.7
  theta <- 1.1
  s <- 0.8

  lambda2 <- 18
  prob2 <- 0.66
  theta2 <- 1.4
  s <- 1.1

  nc <- nimbleCode({
    x ~ dNmixtureAD_BBNB_oneObs(lambda, prob, theta = theta, s = s,
                                Nmin = 0, Nmax = 100)
    prob ~ dunif(0, 1)
    lambda ~ dunif(0, 100)
  })

  Rmodel <- nimbleModel(nc, data = list(x = x),
                        inits = list(prob = prob, theta = theta, s=s,
                                     lambda = lambda),
                        buildDerivs=TRUE)
  Rmodel$calculate()

  Cmodel <- compileNimble(Rmodel)
  Cmodel$calculate()

  nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c('prob', 'lambda', 'theta', 's'))
  v1_case1 <- list(arg1 = c(prob, lambda, theta, s2)) # taping values for prob and lambda
  v2_case1 <- list(arg1 = c(prob2, lambda2, theta2, s2)) # testing values for prob and lambda

  res <- model_calculate_test_case(Rmodel, Cmodel,
                            model_calculate_test, nodesList_case1,
                            v1_case1, v2_case1,
                            0:2, RCrelTol = c(2e-15, 1e-8, 1e-3, 1e-14))
})


test_that("dCJS works with AD", {
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

  res <- model_calculate_test_case(Rmodel, Cmodel,
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

  res <- model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
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

  res <- model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                            nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                            order = 0:2)

######################
#### dCJS_vv case ####

  x <- c(1, 0, 1, 0, 0, 0)
  probSurvive <- c(0.8, 0.5, 0.3, 0.9, 0.9)
  probCapture <- c(1, 0.5, 0.5, 0.4, 0.3, 0.4)

  probSurvive2 <- c(0.7, 0.55, 0.32, 0.8, 0.1)
  probCapture2 <- c(-10, 0.6, 0.7, 0.4, 0.2, 0.2)


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
                                        # v1_case1 <- list(arg1 = c(probSurvive, probCapture)) # taping values for prob and lambda
                                        # v2_case1 <- list(arg1 = c(probSurvive2, probCapture2)) # testing values for prob and lambda

  res <- model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                            nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                            order = 0:2)
})

test_that("dHMM works with AD", {
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

  stuff <- model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                            nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                            order = 0:2)
  #######

  # {
  # ######################
  # # #### dHMM with 0s in transition matrix case ####
  #
  # x <- c(1, 1, 1, 2, 2)
  #
  # init <- c(0.4, 0.2, 0.4)
  # probObs <- t(array(
  #          c(0.9, 0.1,
  #            0.1, 0.9,
  #            0.8, 0.2),
  #          c(2, 3)))
  #
  # probTrans <- t(array(
  #           c(0.3, 0.4, 0.2,
  #             0.1, 0.1, 0.8,
  #             0.05, 0.05, 0.9),
  #           c(3,3)))
  #
  # init2 <- c(0.6, 0.1, 0.3)
  # probObs2 <- t(array(
  #          c(1, 0,
  #            0, 1,
  #            0.7, 0.3),
  #          c(2, 3)))
  #
  # probTrans2 <- t(array(
  #           c(0.4, 0.4, 0.2,
  #             0, 0.3, 0.7,
  #             0, 0, 1),
  #           c(3,3)))
  #
  # nc <- nimbleCode({
  #   x[1:5] ~ dHMM(init[1:3], probObs = probObs[1:3,1:2],
  #                   probTrans = probTrans[1:3, 1:3], len = 5, checkRowSums = 0)
  #   for (i in 1:2) {
  #     init[i] ~ dunif(0, 1)
  #   }
  #   init[3] <- 1 - init[1] - init[2]
  #
  #   for (i in 1:3) {
  #     probObs[i, 1] ~ dunif(0, 1)
  #     probObs[i, 2] <- 1 - probObs[i, 1]
  #     probTrans[i, 1] ~ dunif(0, 1)
  #     probTrans[i, 2] ~ dunif(0, 1 - probTrans[i, 1])
  #     probTrans[i, 3] <- 1 - probTrans[i, 1] - probTrans[i, 2]
  #   }
  #
  # })
  # Rmodel <- nimbleModel(nc, data = list(x = x),
  #                  inits = list(
  #                    init = init,
  #                    probObs = probObs,
  #                    probTrans = probTrans
  #                  ),
  #                  buildDerivs=TRUE)
  # Rmodel$calculate()
  #
  # Cmodel <- compileNimble(Rmodel)
  # Cmodel$calculate()
  #
  # nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c(Rmodel$expandNodeNames('init[1:3]'),
  #                                                                        Rmodel$expandNodeNames('probObs[1:3, 1:2]'),
  #                                                                        Rmodel$expandNodeNames('probTrans[1:3, 1:3]')
  #                                                                        ))
  # v1_case1 <- list(arg1 = c(init[1:3],  probObs[1:3, 1:2], probTrans[1:3, 1:3]))
  # v2_case1 <- list(arg1 = c(init2[1:3], probObs2[1:3, 1:2],  probTrans2[1:3, 1:3]))
  #
  # model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
  #                           nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
  #                           order = 0:2)
  # }
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

  # capture <- capture_warning(
  # The warning is due to getParam for nDim > 2
  Rmodel <- suppressWarnings(nimbleModel(nc, data = list(x = x),
                        inits = list(
                          init = init,
                          probObs = probObs,
                          probTrans = probTrans
                        ),
                        buildDerivs=TRUE)
   )

  Rmodel$calculate()

  Cmodel <- compileNimble(Rmodel)
  Cmodel$calculate()

  nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c(Rmodel$expandNodeNames('init[1:3]'),
                                                                         Rmodel$expandNodeNames('probObs[1:3, 1:2, 1:5]'),
                                                                         Rmodel$expandNodeNames('probTrans[1:3, 1:3]')
                                                                         ))
  v1_case1 <- list(arg1 = c(init[1:3],  probObs[1:3, 1:2, 1:5], probTrans[1:3, 1:3]))
  v2_case1 <- list(arg1 = c(init2[1:3], probObs2[1:3, 1:2, 1:5],  probTrans2[1:3, 1:3]))

  res <- model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                                   nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                                   order = 0:2)
  #######
})

test_that("dDHMM works with AD", {
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

  # capture <- capture_warning(
  Rmodel <- suppressWarnings(nimbleModel(nc, data = list(x = x),
                        inits = list(
                          init = init,
                          probObs = probObs,
                          probTrans = probTrans
                        ),
                        buildDerivs=TRUE)
  )
  Rmodel$calculate()

  Cmodel <- compileNimble(Rmodel)
  Cmodel$calculate()

  nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c(Rmodel$expandNodeNames('init[1:3]'),
                                                                         Rmodel$expandNodeNames('probObs[1:3, 1:2]'),
                                                                         Rmodel$expandNodeNames('probTrans[1:3, 1:3, 1:4]')
                                                                         ))
  v1_case1 <- list(arg1 = c(init[1:3],  probObs[1:3, 1:2], probTrans[1:3, 1:3, 1:4]))
  v2_case1 <- list(arg1 = c(init2[1:3], probObs2[1:3, 1:2],  probTrans2[1:3, 1:3, 1:4]))

  res <- model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                            nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                            order = 0:2,
                            RCrelTol = c(3e-15, 1e-8, 2e-3, 1e-14))
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
  # capture <- capture_warning(
  Rmodel <- suppressWarnings(nimbleModel(nc, data = list(x = x),
                        inits = list(
                          init = init,
                          probObs = probObs,
                          probTrans = probTrans
                        ),
                        buildDerivs=TRUE)
  )
  Rmodel$calculate()

  Cmodel <- compileNimble(Rmodel)
  Cmodel$calculate()

  nodesList_case1 <- setup_update_and_constant_nodes_for_tests(Rmodel, c(Rmodel$expandNodeNames('init[1:3]'),
                                                                         Rmodel$expandNodeNames('probObs[1:3, 1:2, 1:5]'),
                                                                         Rmodel$expandNodeNames('probTrans[1:3, 1:3, 1:4]')
                                                                         ))
  v1_case1 <- list(arg1 = c(init[1:3],  probObs[1:3, 1:2, 1:5], probTrans[1:3, 1:3, 1:4]))
  v2_case1 <- list(arg1 = c(init2[1:3], probObs2[1:3, 1:2, 1:5],  probTrans2[1:3, 1:3, 1:4]))

  res <- model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                            nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                            order = 0:2,
                            RCrelTol = c(1e-14, 1e-8, 2e-3, 2e-14))
  #######
})


test_that("dDynOcc works with AD", {
######################
#### dDynOcc_vvm case ####
                                        # ADtestEnv$RCrelTol sets tolerance
                                        # [1] value [2] first order [3] second order
                                        # can also look at ADtestEnv$CCrelTol

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
  probPersist2 <- c(0.5, 0.55, 0.2)
  probColonize2 <- c(0.3, 0.3, 0.6)
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

  res <- model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                            nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                            order = 0:2,
                            RCrelTol = c(ADtestEnv$RCrelTol[1], 1e-6,
                                         0.004, 1e-7))

######################
#### dDynOcc_vsm case ####

  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)

  init <- 0.7
  probPersist <- c(0.2, 0.1, 0.3)
  probColonize <- 0.4
  p <- matrix(rep(c(0.8, 0.7, 0.8, 0.7, 0.9), each = 4), nrow = 4, byrow =TRUE)

  init2 <- 0.9
  probPersist2 <- c(0.4, 0.4, 0.1)
  probColonize2 <- 0.6
  p2 <- matrix(rep(c(0.7, 0.5, 0.3, 0.8, 0.66), each = 4), nrow = 4, byrow =TRUE)


  nc <- nimbleCode({
    x[1:4, 1:5] ~ dDynOcc_vsm(init,
                              probPersist[1:3],
                              probColonize,
                              p[1:4,1:5],
                              start[1:4], end[1:4])

    init ~ dunif(0, 1)
    probColonize ~ dunif(0, 1)
    for (i in 1:3) {
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
                                                                         Rmodel$expandNodeNames('probColonize'),
                                                                         Rmodel$expandNodeNames('probPersist[1:3]')
                                                                       ))
  v1_case1 <- list(arg1 = c(init,
                            p,
                            probColonize,
                            probPersist
                            ))
  v2_case1 <- list(arg1 = c(init2, p2,
                            probColonize2,
                            probPersist2
                            ))

  res <- model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                            nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                            order = 0:2,
                            RCrelTol = c(ADtestEnv$RCrelTol[1], 4e-7, 0.014, 1e-6))
                                        #0.004, 1e-7))


######################
#### dDynOcc_svm case ####

  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)

  init <- 0.7
  probPersist <- 0.4
  probColonize <- c(0.4, 0.2, 0.1)
  p <- matrix(rep(c(0.8, 0.7, 0.8, 0.8, 0.9), each = 4), nrow = 4, byrow =TRUE)

  init2 <- 0.9
  probPersist2 <- 0.6
  probColonize2 <- c(0.4, 0.2, 0.1)
  p2 <- matrix(rep(c(0.7, 0.5, 0.3, 0.8, 0.66), each = 4), nrow = 4, byrow =TRUE)


  nc <- nimbleCode({
    x[1:4, 1:5] ~ dDynOcc_svm(init,
                              probPersist,
                              probColonize[1:3],
                              p[1:4,1:5],
                              start[1:4], end[1:4])

    init ~ dunif(0, 1)
    for (i in 1:3) {
      probColonize[i] ~ dunif(0, 1)
    }
    probPersist ~ dunif(0, 1)
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
                                                                         Rmodel$expandNodeNames('probPersist')
                                                                       ))
  v1_case1 <- list(arg1 = c(init, p, probColonize, probPersist))
  v2_case1 <- list(arg1 = c(init, p2, probColonize2, probPersist2))

  res <- model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                            nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                            order = 0:2,
                            RCrelTol = c(ADtestEnv$RCrelTol[1], 2e-7, 1e-3, 2e-6))

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

  res <- model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                            nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                            order = 0:2,
                            RCrelTol = c(ADtestEnv$RCrelTol[1], 1e-7, 1e-3, 1e-14))


######################
#### dDynOcc_vsv case ####

  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)

  init <- 0.7
  probPersist <- c(0.4, 0.4, 0.1)
  probColonize <- 0.5
  p <- c(0.8, 0.7, 0.8, 0.8)

  init2 <- 0.9
  probPersist2 <- c(0.4, 0.4, 0.1)
  probColonize2 <- 0.7
  p2 <- c(0.7, 0.5, 0.3, 0.8)


  nc <- nimbleCode({
    x[1:4, 1:5] ~ dDynOcc_vsv(init,
                              probPersist[1:3],
                              probColonize,
                              p[1:4],
                              start[1:4], end[1:4])

    init ~ dunif(0, 1)
    probColonize ~ dunif(0, 1)
    for (i in 1:3) {
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
                                                                         Rmodel$expandNodeNames('probColonize'),
                                                                         Rmodel$expandNodeNames('probPersist[1:3]')
                                                                       ))
  v1_case1 <- list(arg1 = c(init, p, probColonize, probPersist))
  v2_case1 <- list(arg1 = c(init, p2, probColonize2, probPersist2))

  res <- model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                            nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                            order = 0:2,
                            RCrelTol = c(ADtestEnv$RCrelTol[1], 1e-7,
                                         ADtestEnv$RCrelTol[3], 1e-14))
                                        # 0.014, 1e-6))


######################
#### dDynOcc_svv case ####

  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)

  init <- 0.7
  probPersist <- 0.5
  probColonize <- c(0.4, 0.2, 0.1)
  p <- c(0.8, 0.7, 0.8, 0.8)

  init2 <- 0.9
  probPersist2 <- 0.5
  probColonize2 <- c(0.4, 0.2, 0.1)
  p2 <- c(0.7, 0.5, 0.3, 0.8)


  nc <- nimbleCode({
    x[1:4, 1:5] ~ dDynOcc_svv(init,
                              probPersist,
                              probColonize[1:3],
                              p[1:4],
                              start[1:4], end[1:4])

    init ~ dunif(0, 1)
    for (i in 1:3) {
      probColonize[i] ~ dunif(0, 1)
    }
    probPersist ~ dunif(0, 1)
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
                                                                         Rmodel$expandNodeNames('probPersist')
                                                                       ))
  v1_case1 <- list(arg1 = c(init, p, probColonize, probPersist))
  v2_case1 <- list(arg1 = c(init, p2, probColonize2, probPersist2))

  res <- model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                            nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                            order = 0:2)


######################
#### dDynOcc_ssv case ####

  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)

  init <- 0.7
  probPersist <- 0.5
  probColonize <- 0.4
  p <- c(0.8, 0.7, 0.8, 0.8)

  init2 <- 0.9
  probPersist2 <- 0.5
  probColonize2 <- 0.8
  p2 <- c(0.7, 0.5, 0.3, 0.8)


  nc <- nimbleCode({
    x[1:4, 1:5] ~ dDynOcc_ssv(init,
                              probPersist,
                              probColonize,
                              p[1:4],
                              start[1:4], end[1:4])

    init ~ dunif(0, 1)
    probColonize ~ dunif(0, 1)
    probPersist ~ dunif(0, 1)
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
                                                                         Rmodel$expandNodeNames('probColonize'),
                                                                         Rmodel$expandNodeNames('probPersist')
                                                                       ))
  v1_case1 <- list(arg1 = c(init, p, probColonize, probPersist))
  v2_case1 <- list(arg1 = c(init, p2, probColonize2, probPersist2))

  res <- model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
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

  res <- model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                            nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                            order = 0:2)


######################
#### dDynOcc_vss case ####

  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)

  init <- 0.7
  probPersist <- c(0.4, 0.4, 0.1)
  probColonize <- 0.2
  p <- c(0.8)

  init2 <- 0.9
  probPersist2 <- c(0.4, 0.4, 0.1)
  probColonize2 <- 0.8
  p2 <- c(0.7)


  nc <- nimbleCode({
    x[1:4, 1:5] ~ dDynOcc_vss(init,
                              probPersist[1:3],
                              probColonize,
                              p,
                              start[1:4], end[1:4])

    init ~ dunif(0, 1)
    probColonize ~ dunif(0, 1)
    for (i in 1:3) {
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
                                                                         Rmodel$expandNodeNames('probColonize'),
                                                                         Rmodel$expandNodeNames('probPersist[1:3]')
                                                                       ))
  v1_case1 <- list(arg1 = c(init, p, probColonize, probPersist))
  v2_case1 <- list(arg1 = c(init, p2, probColonize2, probPersist2))

  res <- model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                            nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                            order = 0:2,
                            RCrelTol = c(ADtestEnv$RCrelTol[1], ADtestEnv$RCrelTol[2],
                                         ADtestEnv$RCrelTol[3], 2e-14))


######################
#### dDynOcc_svs case ####

  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)

  init <- 0.7
  probPersist <- 0.3
  probColonize <- c(0.4, 0.2, 0.1)
  p <- c(0.8)

  init2 <- 0.9
  probPersist2 <- 0.44
  probColonize2 <- c(0.4, 0.2, 0.1)
  p2 <- c(0.7)


  nc <- nimbleCode({
    x[1:4, 1:5] ~ dDynOcc_svs(init,
                              probPersist,
                              probColonize[1:3],
                              p,
                              start[1:4], end[1:4])

    init ~ dunif(0, 1)
    for (i in 1:3) {
      probColonize[i] ~ dunif(0, 1)
    }
    probPersist ~ dunif(0, 1)
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
                                                                         Rmodel$expandNodeNames('probPersist')
                                                                       ))
  v1_case1 <- list(arg1 = c(init, p, probColonize, probPersist))
  v2_case1 <- list(arg1 = c(init, p2, probColonize2, probPersist2))

  res <- model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                            nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                            order = 0:2)

######################
#### dDynOcc_sss case ####

  x <- matrix(c(0,0,NA,0,
                1,1,1,0,
                0,0,0,0,
                0,0,1,0,
                0,0,0,NA), nrow = 4)
  start <- c(1,1,2,1)
  end <- c(5,5,5,4)

  init <- 0.7
  probPersist <- 0.3
  probColonize <- 0.6
  p <- c(0.8)

  init2 <- 0.9
  probPersist2 <- 0.44
  probColonize2 <- 0.7
  p2 <- c(0.7)


  nc <- nimbleCode({
    x[1:4, 1:5] ~ dDynOcc_sss(init,
                              probPersist,
                              probColonize,
                              p,
                              start[1:4], end[1:4])

    init ~ dunif(0, 1)
    probColonize ~ dunif(0, 1)
    probPersist ~ dunif(0, 1)
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
                                                                         Rmodel$expandNodeNames('probColonize'),
                                                                         Rmodel$expandNodeNames('probPersist')
                                                                       ))
  v1_case1 <- list(arg1 = c(init, p, probColonize, probPersist))
  v2_case1 <- list(arg1 = c(init, p2, probColonize2, probPersist2))

  res <- model_calculate_test_case(Rmodel, Cmodel, deriv_nf = model_calculate_test,
                            nodesList = nodesList_case1, v1 = v1_case1, v2 = v2_case1,
                            order = 0:2)
})

# reset options before finishing
nimbleOptions(enableDerivs = EDopt)
nimbleOptions(buildModelDerivs = BMDopt)
