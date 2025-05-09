skip_on_cran()
skip_if(!requireNamespace("nimbleMacros", quietly=TRUE), 
        "nimbleMacros package unavailable")

set.seed(123)

# Simulate single-species data
M <- 50
J <- 5

x <- rnorm(M)

beta_0 <- 0
beta_x <- 1
alpha_0 <- 0
truth <- c(beta_0, beta_x, alpha_0)

y <- matrix(NA, M, J)

psi <- plogis(beta_0 + beta_x * x)
z <- rbinom(M, 1, psi)
p <- plogis(alpha_0)

for (m in 1:M){
  y[m,] <- rbinom(J, 1, p * z[m])
}

# Simulate multispecies occupancy data
set.seed(123)
S <- 10

x <- rnorm(M)

mu_0 <- 0
sd_0 <- 0.4
beta0 <- rnorm(S, mu_0, sd_0)

mu_x <- 1
sd_x <- 0.3
beta_x <- rnorm(S, mu_x, sd_x)

mu_a <- 0
sd_a <- 0.2
alpha0 <- rnorm(S, mu_a, sd_a)

truth_multi <- c(mu_0, sd_0, mu_x, sd_x, mu_a, sd_a)

ymulti <- array(NA, c(M,J,S))

for (s in 1:S){
  psi <- plogis(beta0[s] + beta_x[s] * x)
  z <- rbinom(M, 1, psi)

  p <- plogis(alpha0[s])

  for (m in 1:M){
    ymulti[m,,s] <- rbinom(J, 1, p * z[m])
  }
}

# To force option to be applied
expect_message(
mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = y, 
                  siteCovs = list(x = x), 
                  obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                  returnModel = TRUE))
mc <- nimbleOptions("enableMacroComments")
nimbleOptions(enableMacroComments = FALSE)
verb <- nimbleOptions("verbose")

test_that("single species occupancy code", {

  mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = y, 
                    siteCovs = list(x = x), 
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    returnModel = TRUE)

  expect_equal(
    mod$getCode(),
    quote({
    for (i_1 in 1:M) {
        logit(psi[i_1]) <- state_Intercept + state_x * x[i_1]
    }
    state_Intercept ~ dunif(-10, 10)
    state_x ~ dlogis(0, 1)
    for (i_2 in 1:M) {
        for (i_3 in 1:J) {
            logit(p[i_2, i_3]) <- det_Intercept
        }
    }
    det_Intercept ~ dunif(-10, 10)
    for (i_4 in 1:M) {
        z[i_4] ~ dbern(psi[i_4])
    }
    for (i_5 in 1:M) {
        for (i_6 in 1:J) {
            y[i_5, i_6] ~ dbern(p[i_5, i_6] * z[i_5])
        }
    }
    })
  )

  pars <- unique(unlist(mod$getMacroParameters()))
  expect_equal(pars, c("psi", "p", "z", "state_Intercept", "state_x",
                       "det_Intercept"))
})

test_that("site covariate in detection formula", {

  mod <- nimbleOccu(stateformula = ~x, detformula = ~x, y = y, 
                    siteCovs = list(x = x), 
                    returnModel = TRUE)

  expect_equal(
    mod$getCode(),
    quote({
    for (i_1 in 1:M) {
        logit(psi[i_1]) <- state_Intercept + state_x * x[i_1]
    }
    state_Intercept ~ dunif(-10, 10)
    state_x ~ dlogis(0, 1)
    for (i_2 in 1:M) {
        for (i_3 in 1:J) {
            logit(p[i_2, i_3]) <- det_Intercept + det_x * x[i_2]
        }
    }
    det_Intercept ~ dunif(-10, 10)
    det_x ~ dlogis(0, 1)
    for (i_4 in 1:M) {
        z[i_4] ~ dbern(psi[i_4])
    }
    for (i_5 in 1:M) {
        for (i_6 in 1:J) {
            y[i_5, i_6] ~ dbern(p[i_5, i_6] * z[i_5])
        }
    }
    })
  )
})

test_that("species cov not allowed in single season model", {
  expect_error(
    nimbleOccu(stateformula = ~x+x3, detformula = ~1, y = y, 
              siteCovs = list(x = x),
              speciesCovs = list(x3 = rnorm(S)),
              returnModel = TRUE),
    "species covariates")
})

test_that("setting custom priors works", {

  pr <- nimbleMacros::setPriors(intercept="dnorm(0, 1)",
                                coefficient = "dnorm(0, 2)")
  dpr <- nimbleMacros::setPriors(intercept="dnorm(0, 3)",
                                coefficient = "dnorm(0, 4)")
  mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = y, 
                    siteCovs = list(x = x), statePriors = pr,
                    detPriors = dpr,
                    returnModel = TRUE)
  expect_equal(
    mod$getCode(),
    quote({
    for (i_1 in 1:M) {
        logit(psi[i_1]) <- state_Intercept + state_x * x[i_1]
    }
    state_Intercept ~ dnorm(0, 1)
    state_x ~ dnorm(0, 2)
    for (i_2 in 1:M) {
        for (i_3 in 1:J) {
            logit(p[i_2, i_3]) <- det_Intercept
        }
    }
    det_Intercept ~ dnorm(0, 3)
    for (i_4 in 1:M) {
        z[i_4] ~ dbern(psi[i_4])
    }
    for (i_5 in 1:M) {
        for (i_6 in 1:J) {
            y[i_5, i_6] ~ dbern(p[i_5, i_6] * z[i_5])
        }
    }
  })
  )
})

test_that("data inputs with wrong dimensions caught", {
  
  x_wrong <- x[-1]

  expect_error(nimbleOccu(stateformula = ~x, detformula = ~1, y = y, 
                    siteCovs = list(x = x_wrong), 
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    returnModel = TRUE),
               "incorrect dimensions")

  x2_wrong <- matrix(rnorm(M*(J-1)), M, J-1)
 
  expect_error(nimbleOccu(stateformula = ~x, detformula = ~x2, y = y, 
                    siteCovs = list(x = x), 
                    obsCovs = list(x2 = x2_wrong),
                    returnModel = TRUE),
               "incorrect dimensions")
})

test_that("marginalized single-species model", {

  mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = y, 
                    siteCovs = list(x = x), 
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    returnModel = TRUE, marginalized = TRUE)

  expect_equal(
    mod$getCode(),
    quote({
    for (i_2 in 1:M) {
        logit(psi[i_2]) <- state_Intercept + state_x * x[i_2]
    }
    state_Intercept ~ dunif(-10, 10)
    state_x ~ dlogis(0, 1)
    for (i_3 in 1:M) {
        for (i_4 in 1:J) {
            logit(p[i_3, i_4]) <- det_Intercept
        }
    }
    det_Intercept ~ dunif(-10, 10)
    for (i_1 in 1:M) {
        y[i_1, 1:J] ~ dOcc_v(probOcc = psi[i_1], probDetect = p[i_1, 1:J], len = J)
    }
    })
  )

  pars <- unique(unlist(mod$getMacroParameters()))
  expect_equal(pars, c("psi", "p", "state_Intercept", "state_x",
                       "det_Intercept"))
})

test_that("silencing model messages", {
  out <- expect_message(mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = y, 
                    siteCovs = list(x = x), 
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    returnModel = TRUE))
  nimbleOptions(verbose=FALSE)
  out <- expect_no_message(mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = y, 
                    siteCovs = list(x = x), 
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    returnModel = TRUE))
  nimbleOptions(verbose=verb)
})

test_that("configurations for single-species latent-state model", {
 
  mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = y, 
                    siteCovs = list(x = x), 
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    returnModel = TRUE)

  # Default RW
  conf <- nimbleEcology:::configureOccuMCMC(mod, "default", samplerControl=list(),
                                            S=1, marginalized=FALSE, occ_var="z") 
  samp <- lapply(conf$getSamplers(), function(x) x$toStr())

  expect_equal(samp,
    list("RW sampler: state_Intercept", "RW sampler: state_x", 
         "RW sampler: det_Intercept", 
         "jointBinary sampler: z[1], z[2], z[4], z[6], z[7], z[8], z[9], z[10], z[11], z[14], z[20], z[21], z[22], z[23], z[24], z[25], z[29], z[31], z[34], z[37], z[40], z[41], z[42], z[43], z[45], z[46], z[47], z[48],  order: TRUE")
  )

  # RW block
  conf <- expect_message(nimbleEcology:::configureOccuMCMC(mod, "RW_block", samplerControl=list(),
                         S=1, marginalized=FALSE, occ_var="z"))
  samp <- lapply(conf$getSamplers(), function(x) x$toStr())

  expect_equal(samp,
    list("RW_block sampler: state_Intercept, state_x, det_Intercept",
         "jointBinary sampler: z[1], z[2], z[4], z[6], z[7], z[8], z[9], z[10], z[11], z[14], z[20], z[21], z[22], z[23], z[24], z[25], z[29], z[31], z[34], z[37], z[40], z[41], z[42], z[43], z[45], z[46], z[47], z[48],  order: TRUE")
  )

  # polya-gamma 
  conf <- nimbleEcology:::configureOccuMCMC(mod, "polyagamma", samplerControl=list(),
                               S=1, marginalized=FALSE, occ_var="z")
  samp <- lapply(conf$getSamplers(), function(x) x$toStr())
  expect_equal(samp,
    list("jointBinary sampler: z[1], z[2], z[4], z[6], z[7], z[8], z[9], z[10], z[11], z[14], z[20], z[21], z[22], z[23], z[24], z[25], z[29], z[31], z[34], z[37], z[40], z[41], z[42], z[43], z[45], z[46], z[47], z[48],  order: TRUE",
    "polyagamma sampler: state_Intercept, state_x,  fixedDesignColumns: TRUE",
    "polyagamma sampler: det_Intercept,  fixedDesignColumns: TRUE"
    )
  )

  # Barker
  conf <- nimbleEcology:::configureOccuMCMC(mod, "barker", samplerControl=list(),
                               S=1, marginalized=FALSE, occ_var="z")
  samp <- lapply(conf$getSamplers(), function(x) x$toStr())
  expect_equal(samp,
    list("barker sampler: state_Intercept, state_x", 
         "barker sampler: det_Intercept", 
         "jointBinary sampler: z[1], z[2], z[4], z[6], z[7], z[8], z[9], z[10], z[11], z[14], z[20], z[21], z[22], z[23], z[24], z[25], z[29], z[31], z[34], z[37], z[40], z[41], z[42], z[43], z[45], z[46], z[47], z[48],  order: TRUE")
  )
})

test_that("mcmc configurations for single-species marginalized model", {

  mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = y, 
                    siteCovs = list(x = x), 
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    returnModel = TRUE, marginalized = TRUE, buildDerivs=TRUE)
  # RW
  conf <- nimbleEcology:::configureOccuMCMC(mod, "default", samplerControl=list(),
                               S=1, marginalized=TRUE) 
  samp <- lapply(conf$getSamplers(), function(x) x$toStr())
  expect_equal(samp,
    list("RW sampler: state_Intercept", "RW sampler: state_x", 
         "RW sampler: det_Intercept")
  )

  # RW_block
  conf <- expect_message(nimbleEcology:::configureOccuMCMC(mod, "RW_block", 
                          samplerControl=list(), S=1, marginalized=TRUE))
  samp <- lapply(conf$getSamplers(), function(x) x$toStr())
  expect_equal(samp,
    list("RW_block sampler: state_Intercept, state_x, det_Intercept")
  )
  
  # barker
  conf <- nimbleEcology:::configureOccuMCMC(mod, "barker", samplerControl=list(),
                               S=1, marginalized=TRUE)
  samp <- lapply(conf$getSamplers(), function(x) x$toStr())
  expect_equal(samp,
    list("barker sampler: state_Intercept, state_x", "barker sampler: det_Intercept")
  )

  # HMC
  conf <- nimbleEcology:::configureOccuMCMC(mod, "hmc", samplerControl=list(),
                               S=S, marginalized=TRUE)
  samp <- lapply(conf$getSamplers(), function(x) x$toStr())
  expect_equal(samp,
    list("NUTS sampler: state_Intercept, state_x, det_Intercept,  warmupMode: default")
  )
})

test_that("fit single-species latent state model", {
  nimbleOptions(verbose=FALSE)
  set.seed(123)
  mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = y, 
                    siteCovs = list(x = x), 
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    nburnin=1000, niter=2000, nchains=2,
                    samplesAsCodaMCMC = TRUE)

  mat <- as.matrix(mod)
  expect_equal(ncol(mat), M+3) 
  mat <- mat[,grepl("_", colnames(mat), fixed=TRUE)]
  est <- round(colMeans(mat), 2)
  expect_equivalent(est, c(-0.14, -0.18, 0.85))
  
  # Also save p
  mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = y, 
                    siteCovs = list(x = x), 
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    nburnin=1000, niter=2000, nchains=2, saveP=TRUE,
                    samplesAsCodaMCMC = TRUE)
  mat <- as.matrix(mod)
  expect_equal(ncol(mat), 3+M+J*M)
  nimbleOptions(verbose=verb)
})

test_that("fit single-species latent state model with polyagamma", {
  nimbleOptions(verbose=FALSE)
  set.seed(123)
  pr <- nimbleMacros::setPriors(intercept="dnorm(0, sd=10)",
                                coefficient="dnorm(0, sd=10)")
  mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = y, 
                    siteCovs = list(x = x), 
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    statePriors = pr, detPriors = pr,
                    nburnin=1000, niter=2000, nchains=2,
                    sampler = "polyagamma",
                    samplesAsCodaMCMC = TRUE)

  mat <- as.matrix(mod)
  expect_equal(ncol(mat), M+3) 
  mat <- mat[,grepl("_", colnames(mat), fixed=TRUE)]
  est <- round(colMeans(mat), 2)
  expect_equivalent(est, c(-0.14, -0.19, 0.88))
  nimbleOptions(verbose=verb)
})

test_that("fit single-species latent state model with barker", {
  nimbleOptions(verbose=FALSE)
  set.seed(123)
  mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = y, 
                    siteCovs = list(x = x), 
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    nburnin=1000, niter=2000, nchains=2,
                    sampler = "barker",
                    samplesAsCodaMCMC = TRUE)

  mat <- as.matrix(mod)
  expect_equal(ncol(mat), M+3) 
  mat <- mat[,grepl("_", colnames(mat), fixed=TRUE)]
  est <- round(colMeans(mat), 2)
  expect_equivalent(est, c(-0.13, -0.18, 0.82))
  nimbleOptions(verbose=verb)
})

test_that("fit single-species marginalized model with hmc", {
  nimbleOptions(verbose=FALSE)
  set.seed(123)
  mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = y, 
                    siteCovs = list(x = x), 
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    nburnin=1000, niter=2000, nchains=2,
                    samplesAsCodaMCMC = TRUE, marginalized=TRUE, sampler="hmc")

  mat <- as.matrix(mod)
  expect_equal(ncol(mat), M+3) 
  mat <- mat[,grepl("_", colnames(mat), fixed=TRUE)]
  est <- round(colMeans(mat), 2)
  expect_equivalent(est, c(-0.13, -0.19, 0.83))
  nimbleOptions(verbose=verb)
})

nimbleOptions(enableMacroComments = mc)
nimbleOptions(verbose = verb)
