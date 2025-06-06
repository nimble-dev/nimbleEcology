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
                               S=1, marginalized=TRUE)
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

test_that("multi species occupancy code", {

  mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = ymulti, 
                    siteCovs = list(x = x), 
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    returnModel = TRUE)

  expect_equal(
    mod$getCode(),
    quote({
    for (i_1 in 1:M) {
        for (i_2 in 1:S) {
            logit(psi[i_1, i_2]) <- state_speciesID[speciesID[i_2]] + 
                state_x_speciesID[speciesID[i_2]] * x[i_1]
        }
    }
    state_Intercept ~ dunif(-10, 10)
    state_x ~ dlogis(0, 1)
    state_sd_speciesID ~ dhalfflat()
    for (i_3 in 1:10) {
        state_speciesID[i_3] ~ dnorm(state_Intercept, sd = state_sd_speciesID)
    }
    state_sd_x_speciesID ~ dhalfflat()
    for (i_4 in 1:10) {
        state_x_speciesID[i_4] ~ dnorm(state_x, sd = state_sd_x_speciesID)
    }
    for (i_5 in 1:M) {
        for (i_6 in 1:J) {
            for (i_7 in 1:S) {
                logit(p[i_5, i_6, i_7]) <- det_speciesID[speciesID[i_7]]
            }
        }
    }
    det_Intercept ~ dunif(-10, 10)
    det_sd_speciesID ~ dhalfflat()
    for (i_8 in 1:10) {
        det_speciesID[i_8] ~ dnorm(det_Intercept, sd = det_sd_speciesID)
    }
    for (i_9 in 1:M) {
        for (i_10 in 1:S) {
            z[i_9, i_10] ~ dbern(psi[i_9, i_10])
        }
    }
    for (i_11 in 1:M) {
        for (i_12 in 1:J) {
            for (i_13 in 1:S) {
                y[i_11, i_12, i_13] ~ dbern(p[i_11, i_12, i_13] * 
                  z[i_11, i_13])
            }
        }
    }
    })
  )

  pars <- unique(unlist(mod$getMacroParameters()))
  expect_equal(pars, c("psi", "p", "z", "speciesID", "state_speciesID", "state_x_speciesID", 
    "det_speciesID", "state_Intercept", "state_x", "state_sd_speciesID", 
    "state_sd_x_speciesID", "det_Intercept", "det_sd_speciesID"))
})

test_that("species cov in model", {
  mod <- nimbleOccu(stateformula = ~x+x3, detformula = ~1, y = ymulti, 
              siteCovs = list(x = x),
              speciesCovs = list(x3 = rnorm(S)),
              returnModel = TRUE)

  expect_equal(
    mod$getCode(),
    quote({
    for (i_1 in 1:M) {
        for (i_2 in 1:S) {
            logit(psi[i_1, i_2]) <- state_x3 * x3[i_2] + state_speciesID[speciesID[i_2]] + 
                state_x_speciesID[speciesID[i_2]] * x[i_1]
        }
    }
    state_Intercept ~ dunif(-10, 10)
    state_x ~ dlogis(0, 1)
    state_x3 ~ dlogis(0, 1)
    state_sd_speciesID ~ dhalfflat()
    for (i_3 in 1:10) {
        state_speciesID[i_3] ~ dnorm(state_Intercept, sd = state_sd_speciesID)
    }
    state_sd_x_speciesID ~ dhalfflat()
    for (i_4 in 1:10) {
        state_x_speciesID[i_4] ~ dnorm(state_x, sd = state_sd_x_speciesID)
    }
    for (i_5 in 1:M) {
        for (i_6 in 1:J) {
            for (i_7 in 1:S) {
                logit(p[i_5, i_6, i_7]) <- det_speciesID[speciesID[i_7]]
            }
        }
    }
    det_Intercept ~ dunif(-10, 10)
    det_sd_speciesID ~ dhalfflat()
    for (i_8 in 1:10) {
        det_speciesID[i_8] ~ dnorm(det_Intercept, sd = det_sd_speciesID)
    }
    for (i_9 in 1:M) {
        for (i_10 in 1:S) {
            z[i_9, i_10] ~ dbern(psi[i_9, i_10])
        }
    }
    for (i_11 in 1:M) {
        for (i_12 in 1:J) {
            for (i_13 in 1:S) {
                y[i_11, i_12, i_13] ~ dbern(p[i_11, i_12, i_13] * 
                  z[i_11, i_13])
            }
        }
    }
  })
  )

  pars <- unique(unlist(mod$getMacroParameters()))
  expect_equal(pars, c("psi", "p", "z", "speciesID", "state_x3", "state_speciesID", "state_x_speciesID", 
    "det_speciesID", "state_Intercept", "state_x", "state_sd_speciesID", 
    "state_sd_x_speciesID", "det_Intercept", "det_sd_speciesID"))
})

test_that("error when species cov is incorrect length", {
  x3_wrong <- rnorm(S-1)
  expect_error(nimbleOccu(stateformula = ~x+x3, detformula = ~1, y = ymulti, 
              siteCovs = list(x = x),
              speciesCovs = list(x3 = x3_wrong),
              returnModel = TRUE), "incorrect dimensions")
})

test_that("marginalized multi-species model", {

  mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = ymulti, 
                    siteCovs = list(x = x), 
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    returnModel = TRUE, marginalized = TRUE)

  expect_equal(
    mod$getCode(),
    quote({
    for (i_3 in 1:M) {
        for (i_4 in 1:S) {
            logit(psi[i_3, i_4]) <- state_speciesID[speciesID[i_4]] + 
                state_x_speciesID[speciesID[i_4]] * x[i_3]
        }
    }
    state_Intercept ~ dunif(-10, 10)
    state_x ~ dlogis(0, 1)
    state_sd_speciesID ~ dhalfflat()
    for (i_5 in 1:10) {
        state_speciesID[i_5] ~ dnorm(state_Intercept, sd = state_sd_speciesID)
    }
    state_sd_x_speciesID ~ dhalfflat()
    for (i_6 in 1:10) {
        state_x_speciesID[i_6] ~ dnorm(state_x, sd = state_sd_x_speciesID)
    }
    for (i_7 in 1:M) {
        for (i_8 in 1:J) {
            for (i_9 in 1:S) {
                logit(p[i_7, i_8, i_9]) <- det_speciesID[speciesID[i_9]]
            }
        }
    }
    det_Intercept ~ dunif(-10, 10)
    det_sd_speciesID ~ dhalfflat()
    for (i_10 in 1:10) {
        det_speciesID[i_10] ~ dnorm(det_Intercept, sd = det_sd_speciesID)
    }
    for (i_1 in 1:M) {
        for (i_2 in 1:S) {
            y[i_1, 1:J, i_2] ~ dOcc_v(probOcc = psi[i_1, i_2], 
                probDetect = p[i_1, 1:J, i_2], len = J)
        }
    }
    })
  )

  pars <- unique(unlist(mod$getMacroParameters()))
  expect_equal(pars, c("psi", "p", "speciesID", "state_speciesID", "state_x_speciesID", 
    "det_speciesID", "state_Intercept", "state_x", "state_sd_speciesID", 
    "state_sd_x_speciesID", "det_Intercept", "det_sd_speciesID"))
})

test_that("mcmc configurations for multi-species latent state model", {
  mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = ymulti, 
                    siteCovs = list(x = x), 
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    returnModel = TRUE)

  # Default RW
  conf <- nimbleEcology:::configureOccuMCMC(mod, "default", samplerControl=list(),
                                            S=S, marginalized=FALSE, occ_var="z") 
  samp <- lapply(conf$getSamplers(), function(x) x$toStr())
  expect_equal(length(samp), 37)
  expect_equal(samp[[1]], "RW sampler: state_Intercept")
  expect_equal(samp[[3]], "conjugate_dhalfflat_dnorm_identity sampler: state_sd_speciesID")
  expect_equal(samp[[7]], "RW sampler: state_speciesID[1]")
  expect_equal(samp[[17]], "RW sampler: state_x_speciesID[1]")
  expect_equal(samp[[27]], "RW sampler: det_speciesID[1]")
  expect_true(grepl("jointBinary", samp[[37]]))

  # RW block
  conf <- expect_message(nimbleEcology:::configureOccuMCMC(mod, "RW_block", samplerControl=list(),
                         S=S, marginalized=FALSE, occ_var="z"))
  samp <- lapply(conf$getSamplers(), function(x) x$toStr())
  
  expect_equal(length(samp), 17)
  expect_equal(samp[[1]], "RW sampler: state_Intercept")
  expect_equal(samp[[3]], "conjugate_dhalfflat_dnorm_identity sampler: state_sd_speciesID")
  expect_equal(samp[[7]], "RW_block sampler: state_speciesID[1], state_x_speciesID[1], det_speciesID[1]")
  expect_true(grepl("jointBinary", samp[[17]]))

  # polya-gamma 
  conf <- nimbleEcology:::configureOccuMCMC(mod, "polyagamma", samplerControl=list(),
                               S=S, marginalized=FALSE, occ_var="z")
  samp <- lapply(conf$getSamplers(), function(x) x$toStr())

  expect_equal(length(samp), 27)
  expect_equal(samp[[1]], "RW sampler: state_Intercept")
  expect_equal(samp[[3]], "conjugate_dhalfflat_dnorm_identity sampler: state_sd_speciesID")
  expect_true(grepl("jointBinary", samp[[7]]))
  expect_equal(samp[[8]], "polyagamma sampler: state_speciesID[1], state_x_speciesID[1],  fixedDesignColumns: TRUE") 
  expect_equal(samp[[9]], "polyagamma sampler: det_speciesID[1],  fixedDesignColumns: TRUE")

  # Barker
  conf <- nimbleEcology:::configureOccuMCMC(mod, "barker", samplerControl=list(),
                               S=S, marginalized=FALSE, occ_var="z")
  samp <- lapply(conf$getSamplers(), function(x) x$toStr())

  expect_equal(length(samp), 27)
  expect_equal(samp[[1]], "RW sampler: state_Intercept")
  expect_equal(samp[[3]], "conjugate_dhalfflat_dnorm_identity sampler: state_sd_speciesID")
  expect_equal(samp[[7]], "barker sampler: state_speciesID[1], state_x_speciesID[1]") 
  expect_equal(samp[[8]], "barker sampler: det_speciesID[1]")
  expect_true(grepl("jointBinary", samp[[27]]))
})

test_that("mcmc configurations for multi-species marginalized model", {

  mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = ymulti, 
                    siteCovs = list(x = x), 
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    returnModel = TRUE, marginalized=TRUE, buildDerivs=TRUE)
  # RW
  conf <- nimbleEcology:::configureOccuMCMC(mod, "default", samplerControl=list(),
                               S=S, marginalized=TRUE) 
  samp <- lapply(conf$getSamplers(), function(x) x$toStr())
  expect_equal(length(samp), 36)
  expect_equal(samp[[1]], "RW sampler: state_Intercept")
  expect_equal(samp[[3]], "conjugate_dhalfflat_dnorm_identity sampler: state_sd_speciesID")
  expect_equal(samp[[7]], "RW sampler: state_speciesID[1]")
  expect_equal(samp[[17]], "RW sampler: state_x_speciesID[1]")
  expect_equal(samp[[27]], "RW sampler: det_speciesID[1]")

  # RW_block
  conf <- expect_message(nimbleEcology:::configureOccuMCMC(mod, "RW_block", 
                          samplerControl=list(), S=S, marginalized=TRUE))
  samp <- lapply(conf$getSamplers(), function(x) x$toStr())
  expect_equal(length(samp), 16)
  expect_equal(samp[[1]], "RW sampler: state_Intercept")
  expect_equal(samp[[3]], "conjugate_dhalfflat_dnorm_identity sampler: state_sd_speciesID")
  expect_equal(samp[[7]], "RW_block sampler: state_speciesID[1], state_x_speciesID[1], det_speciesID[1]")
  
  # barker
  conf <- nimbleEcology:::configureOccuMCMC(mod, "barker", samplerControl=list(),
                               S=S, marginalized=TRUE)
  samp <- lapply(conf$getSamplers(), function(x) x$toStr())
  expect_equal(length(samp), 26)
  expect_equal(samp[[1]], "RW sampler: state_Intercept")
  expect_equal(samp[[3]], "conjugate_dhalfflat_dnorm_identity sampler: state_sd_speciesID")
  expect_equal(samp[[7]], "barker sampler: state_speciesID[1], state_x_speciesID[1]") 
  expect_equal(samp[[8]], "barker sampler: det_speciesID[1]")

  # HMC
  conf <- nimbleEcology:::configureOccuMCMC(mod, "hmc", samplerControl=list(),
                               S=S, marginalized=TRUE)
  samp <- lapply(conf$getSamplers(), function(x) x$toStr())
  expect_equal(length(samp), 1)
  expect_true(grepl("NUTS", samp[[1]]))
})

test_that("fit multi-species latent state model", {
  nimbleOptions(verbose=FALSE)
  set.seed(123)
  pr <- nimbleMacros::setPriors(sd="dhalfflat()", intercept="dunif(-10, 10)", coefficient="dnorm(0, sd=10)")
  mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = ymulti, 
                    siteCovs = list(x = x),
                    statePriors=pr, detPriors=pr,
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    nburnin=3000, niter=4000, nchains=1,
                    samplesAsCodaMCMC = TRUE)

  mat <- as.matrix(mod)
  expect_equal(ncol(mat), 506) 
  mat <- mat[,grepl("_", colnames(mat), fixed=TRUE)]
  est <- round(colMeans(mat), 2)
  expect_equivalent(est, c(-0.18, 0.19, 0.16, 0.23, 0.23, 1.06))
  nimbleOptions(verbose=verb)
})

test_that("fit single-species latent state model with polyagamma", {
  nimbleOptions(verbose=FALSE)
  nimbleOptions(verbose=FALSE)
  set.seed(123)
  pr <- nimbleMacros::setPriors(sd="dhalfflat()", intercept="dnorm(0, sd=10)", coefficient="dnorm(0, sd=10)")
  mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = ymulti, 
                    siteCovs = list(x = x),
                    statePriors=pr, detPriors=pr,
                    sampler = "polyagamma",
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    nburnin=3000, niter=4000, nchains=1,
                    samplesAsCodaMCMC = TRUE)

  mat <- as.matrix(mod)
  expect_equal(ncol(mat), 506) 
  mat <- mat[,grepl("_", colnames(mat), fixed=TRUE)]
  est <- round(colMeans(mat), 2)
  expect_equivalent(est, c(-0.18, 0.16, 0.13, 0.16, 0.21, 1.00))
  nimbleOptions(verbose=verb)
})

test_that("fit single-species latent state model with barker", {
  nimbleOptions(verbose=FALSE)
  set.seed(123)
  pr <- nimbleMacros::setPriors(sd="dunif(0, 3)", intercept="dunif(-5, 5)", coefficient="dnorm(0, sd=2.5)")
  mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = ymulti, 
                    siteCovs = list(x = x),
                    statePriors=pr, detPriors=pr,
                    sampler="barker",
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    nburnin=3000, niter=4000, nchains=1,
                    samplesAsCodaMCMC = TRUE)

  mat <- as.matrix(mod)
  expect_equal(ncol(mat), 506) 
  mat <- mat[,grepl("_", colnames(mat), fixed=TRUE)]
  est <- round(colMeans(mat), 2)
  expect_equivalent(est, c(-0.18, 0.20, 0.16, 0.23, 0.24, 1.03))
  nimbleOptions(verbose=verb)
})

test_that("fit single-species marginalized model with hmc", {
  nimbleOptions(verbose=FALSE)
  set.seed(123)
  mod <- nimbleOccu(stateformula = ~x, detformula = ~1, y = ymulti, 
                    siteCovs = list(x = x), 
                    obsCovs = list(x2 = matrix(rnorm(M*J), M, J)),
                    nburnin=300, niter=500, nchains=1,
                    samplesAsCodaMCMC = TRUE, marginalized=TRUE, sampler="hmc")

  mat <- as.matrix(mod)
  expect_equal(ncol(mat), 506) 
  mat <- mat[,grepl("_", colnames(mat), fixed=TRUE)]
  est <- round(colMeans(mat), 2)
  expect_equivalent(est, c(-0.18, 0.19, 0.14, 0.25, 0.20, 1.0))
  nimbleOptions(verbose=verb)
})


nimbleOptions(enableMacroComments = mc)
nimbleOptions(verbose = verb)
