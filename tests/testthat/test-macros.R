skip_on_cran()
nimbleMacros_installed <- requireNamespace("nimbleMacros")
skip_if_not(nimbleMacros_installed)

macro_setting <- nimbleOptions("enableModelMacros")
comm_setting <- nimbleOptions("enableMacroComments")
nimbleOptions(enableModelMacros = TRUE)

# Simulate OCCUPANCY dataset
set.seed(123)

M <- 100
J <- 3
x <- rnorm(M)
b0 <- 0
b1 <- 0.5
psi <- plogis(b0 + b1 * x)
z <- rbinom(M, 1, psi)
g <- sample(letters[1:8], M, replace=TRUE)
x2  <- matrix(rnorm(M*J), M, J)
a0 <- 0
a1 <- -0.5
p <- plogis(a0 + a1*x2)
y <- matrix(NA, M, J)
for (i in 1:M){
  for (j in 1:J){
    y[i,j] <- rbinom(1, 1, p[i,j]*z[i])
  }
}
# Add some NAs
y[1,1] <- NA
y[2,2:3] <- NA
y[3,] <- NA

# the macro really needs to be able to modify data...

# Problem: macro creates z in constants
# then modelDef moves z to data (correctly)
# then nimbleModel does a separate check on the *original* constants,
# which of course now are missing z, and so z doesn't end up in the final
# constants and is all NA

#z <- apply(y, 1, max, na.rm=TRUE)
#z[z==0] <- NA
const <- list(y=y, x=x, x2=x2, M=M, J=J, g = g)

sc <- list(x = x)
oc <- list(x2 = x2)

# Force nimbleMacros to load namespace so we can turn off comments
mod <- nimbleOccu(~1, ~1, y=y, returnModel = TRUE)
nimbleOptions(enableMacroComments = FALSE)

test_that("Latent-state OCCUPANCY macro works", {

  # Latent state, no covariates
  code <- nimbleCode({
    y[1:M, 1:J] ~ OCCUPANCY(~1, ~1)
  })

  mod <- nimbleModel(code, constants=const)

  code_ref <- quote({
    for (i_1 in 1:M) {
        logit(psi[i_1]) <- state_Intercept
    }
    state_Intercept ~ dunif(-10, 10)
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
  
  expect_equal(mod$getCode(), code_ref)

  # Make sure z is properly initialized
  expect_equal(mod$z, apply(const$y, 1, function(x) if(all(is.na(x))) return(0) else max(x, na.rm=TRUE)))

  # With covariates (unbracketed)
  code <- nimbleCode({
    y[1:M, 1:J] ~ OCCUPANCY(~x, ~x2)
  })

  code_ref <- quote({
    for (i_1 in 1:M) {
        logit(psi[i_1]) <- state_Intercept + state_x * x[i_1]
    }
    state_Intercept ~ dunif(-10, 10)
    state_x ~ dlogis(0, 1)
    for (i_2 in 1:M) {
        for (i_3 in 1:J) {
            logit(p[i_2, i_3]) <- det_Intercept + det_x2 * x2[i_2,
                i_3]
        }
    }
    det_Intercept ~ dunif(-10, 10)
    det_x2 ~ dlogis(0, 1)
    for (i_4 in 1:M) {
        z[i_4] ~ dbern(psi[i_4])
    }
    for (i_5 in 1:M) {
        for (i_6 in 1:J) {
            y[i_5, i_6] ~ dbern(p[i_5, i_6] * z[i_5])
        }
    }
  })
  
  mod <- nimbleModel(code, constants=const)
  
  expect_equal(mod$getCode(), code_ref)

  # With covariates (bracketed)
  code <- nimbleCode({
    y[1:M, 1:J] ~ OCCUPANCY(~x[1:M], ~x2[1:M,1:J])
  })
  
  mod <- nimbleModel(code, constants=const)
  
  expect_equal(mod$getCode(), code_ref)

})


test_that("Marginalized OCCUPANCY macro works", {
  # Marginalized version of the model
  code <- nimbleCode({
    y[1:M, 1:J] ~ OCCUPANCY(~x, ~x2, marginalized=TRUE)
  })

  code_ref <- quote({
    for (i_2 in 1:M) {
        logit(psi[i_2]) <- state_Intercept + state_x * x[i_2]
    }
    state_Intercept ~ dunif(-10, 10)
    state_x ~ dlogis(0, 1)
    for (i_3 in 1:M) {
        for (i_4 in 1:J) {
            logit(p[i_3, i_4]) <- det_Intercept + det_x2 * x2[i_3,
                i_4]
        }
    }
    det_Intercept ~ dunif(-10, 10)
    det_x2 ~ dlogis(0, 1)
    for (i_1 in 1:M) {
        y[i_1, 1:J] ~ dOcc_v(probOcc = psi[i_1], probDetect = p[i_1,
            1:J], len = J)
    }
  })

  mod <- nimbleModel(code, constants=const)

  expect_equal(mod$getCode(), code_ref)

})

test_that("partially centered random effect works", {

  code <- nimbleCode({
    y[1:M, 1:J] ~ OCCUPANCY(~x + (1|g), ~x2, centerVar=g)
  })

  code_ref <- quote({
    for (i_1 in 1:M) {
        logit(psi[i_1]) <- state_x * x[i_1] + state_g[g[i_1]]
    }
    state_Intercept ~ dunif(-10, 10)
    state_x ~ dlogis(0, 1)
    state_sd_g ~ dunif(0, 100)
    for (i_2 in 1:8) {
        state_g[i_2] ~ dnorm(state_Intercept, sd = state_sd_g)
    }
    for (i_3 in 1:M) {
        for (i_4 in 1:J) {
            logit(p[i_3, i_4]) <- det_Intercept + det_x2 * x2[i_3,
                i_4]
        }
    }
    det_Intercept ~ dunif(-10, 10)
    det_x2 ~ dlogis(0, 1)
    for (i_5 in 1:M) {
        z[i_5] ~ dbern(psi[i_5])
    }
    for (i_6 in 1:M) {
        for (i_7 in 1:J) {
            y[i_6, i_7] ~ dbern(p[i_6, i_7] * z[i_6])
        }
    }
  })

  mod <- nimbleModel(code, constants=const)

  expect_equal(mod$getCode(), code_ref)
})

test_that("varying sampling occasions works", {

  K <- rep(J, M)
  const <- list(y=y, x=x, x2=x2, M=M, K=K, g = g)
  code <- nimbleCode({
    y[1:M, 1:K[1:M]] ~ OCCUPANCY(~1, ~x2[1:M, 1:K[1:M]])
  })

  code_ref <- quote({
    for (i_1 in 1:M) {
        logit(psi[i_1]) <- state_Intercept
    }
    state_Intercept ~ dunif(-10, 10)
    for (i_2 in 1:M) {
        for (i_3 in 1:K[i_2]) {
            logit(p[i_2, i_3]) <- det_Intercept + det_x2 * x2[i_2,
                i_3]
        }
    }
    det_Intercept ~ dunif(-10, 10)
    det_x2 ~ dlogis(0, 1)
    for (i_4 in 1:M) {
        z[i_4] ~ dbern(psi[i_4])
    }
    for (i_5 in 1:M) {
        for (i_6 in 1:K[i_5]) {
            y[i_5, i_6] ~ dbern(p[i_5, i_6] * z[i_5])
        }
    }
  })


  mod <- nimbleModel(code, constants=const)

  expect_equal(mod$getCode(), code_ref)
})

test_that("multispeciesFormula works with species covs", {

  expect_equal(
    multispeciesFormula(~x[1:M] + z[1:S], quote(1:S), quote(sp)),
    ~x[1:M] + z[1:S] + (1 + x[1:M] + z[1:S] || sp[1:S])
  )

  expect_equal(
    multispeciesFormula(~x[1:M] + z[1:S], quote(1:S), quote(sp), "z"),
    ~x[1:M] + z[1:S] + (1 + x[1:M] || sp[1:S])
  )
  
  expect_equal(
    multispeciesFormula(~x[1:M]*z[1:S], quote(1:S), quote(sp), "z"),
    ~x[1:M]*z[1:S] + (1 + x[1:M] || sp[1:S])
  )
})

test_that("nimbleOccu works", {
  sc <- list(x = x)
  oc <- list(x2 = x2)
  mod <- nimbleOccu(~x, ~x2, y, siteCovs=sc, obsCovs=oc, returnModel=TRUE)

  code_ref <- quote({
    for (i_1 in 1:M) {
        logit(psi[i_1]) <- state_Intercept + state_x * x[i_1]
    }
    state_Intercept ~ dunif(-10, 10)
    state_x ~ dlogis(0, 1)
    for (i_2 in 1:M) {
        for (i_3 in 1:J) {
            logit(p[i_2, i_3]) <- det_Intercept + det_x2 * x2[i_2,
                i_3]
        }
    }
    det_Intercept ~ dunif(-10, 10)
    det_x2 ~ dlogis(0, 1)
    for (i_4 in 1:M) {
        z[i_4] ~ dbern(psi[i_4])
    }
    for (i_5 in 1:M) {
        for (i_6 in 1:J) {
            y[i_5, i_6] ~ dbern(p[i_5, i_6] * z[i_5])
        }
    }
  })
  expect_equal(mod$getCode(), code_ref)

  fit <- nimbleOccu(~x, ~x2, y, siteCovs=sc, obsCovs=oc, nchain=2,
                    niter=100, nburnin=50)
  expect_is(fit, "list")
  expect_equal(dim(fit$chain1), c(50, 104))

  # Different dimensions for different obs covs
  mod <- nimbleOccu(~1, ~x2 + x, y, siteCovs=sc, obsCovs=oc, returnModel=TRUE)
  code_ref <- quote({
    for (i_1 in 1:M) {
        logit(psi[i_1]) <- state_Intercept
    }
    state_Intercept ~ dunif(-10, 10)
    for (i_2 in 1:M) {
        for (i_3 in 1:J) {
            logit(p[i_2, i_3]) <- det_Intercept + det_x2 * x2[i_2,
                i_3] + det_x * x[i_2]
        }
    }
    det_Intercept ~ dunif(-10, 10)
    det_x2 ~ dlogis(0, 1)
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
  expect_equal(code_ref, mod$getCode())

  # Check marginalized model
  mod <- nimbleOccu(~x, ~x2, y, siteCovs=sc, obsCovs=oc, returnModel=TRUE,
                    marginalized=TRUE)
  docc_code <- mod$getCode()[[8]]
  expect_equal(docc_code,
    quote({
      for (i_1 in 1:M) {
        y[i_1, 1:J] ~ dOcc_v(probOcc = psi[i_1], probDetect = p[i_1,
          1:J], len = J)
      }
    })[[2]]
  )

  # Check null model
  mod <- nimbleOccu(~1, ~1, y, siteCovs=sc, obsCovs=oc, returnModel=TRUE)

  # Can't have species covs in a single-species model
  sp <- list(x3 = rnorm(3))
  expect_error(nimbleOccu(~1, ~1, y, siteCovs=sc, obsCovs=oc, speciesCovs=sp, returnModel=TRUE))
})

# MULTISPECIESOCCUPANCY--------------------------------------------------------

# Simulate data
M <- 100
S <- 10
J <- 5

a0_mn <- 0
a0_sd <- 0.2
a1_mn <- -0.3
a1_sd <- 0.1

b0_mn <- 0
b0_sd <- 0.1
b1_mn <- 0.2
b1_sd <- 0.1

set.seed(1)

a0 <- rnorm(S, a0_mn, a0_sd)
a1 <- rnorm(S, a1_mn, a1_sd)
b0 <- rnorm(S, b0_mn, b0_sd)
b1 <- rnorm(S, b1_mn, b1_sd)

x1 <- rnorm(M)
x2 <- matrix(rnorm(M*J), M, J)

y <- array(NA, c(M, J, S))

for (s in 1:S){
  psi <- plogis(a0[s] + a1[s] * x1)
  z <- rbinom(M, 1, psi)
  p <- matrix(NA, M, J)
  for (n in 1:M){
    p[n,] <- plogis(b0[s] + b1[s] * x2[n,])
    y[n,,s] <- rbinom(J, 1, p[n,] * z[n])
  }
}
constants <- list(y = y, M = M, J = J, S = S, x1 = x1, x2 = x2)

code_ref <- quote({
  for (i_1 in 1:M) {
      for (i_2 in 1:S) {
          logit(psi[i_1, i_2]) <- state_speciesID[speciesID[i_2]] +
              state_speciesID_x1[speciesID[i_2]] * x1[i_1]
      }
  }
  state_Intercept ~ dunif(-10, 10)
  state_x1 ~ dlogis(0, 1)
  state_sd_speciesID ~ dunif(0, 5)
  for (i_3 in 1:10) {
      state_speciesID[i_3] ~ dnorm(state_Intercept, sd = state_sd_speciesID)
  }
  state_sd_speciesID_x1 ~ dunif(0, 5)
  for (i_4 in 1:10) {
      state_speciesID_x1[i_4] ~ dnorm(state_x1, sd = state_sd_speciesID_x1)
  }
  for (i_5 in 1:M) {
      for (i_6 in 1:J) {
          for (i_7 in 1:S) {
              logit(p[i_5, i_6, i_7]) <- det_speciesID[speciesID[i_7]] +
                det_speciesID_x2[speciesID[i_7]] * x2[i_5,
                  i_6]
          }
      }
  }
  det_Intercept ~ dunif(-10, 10)
  det_x2 ~ dlogis(0, 1)
  det_sd_speciesID ~ dunif(0, 5)
  for (i_8 in 1:10) {
      det_speciesID[i_8] ~ dnorm(det_Intercept, sd = det_sd_speciesID)
  }
  det_sd_speciesID_x2 ~ dunif(0, 5)
  for (i_9 in 1:10) {
      det_speciesID_x2[i_9] ~ dnorm(det_x2, sd = det_sd_speciesID_x2)
  }
  for (i_10 in 1:M) {
      for (i_11 in 1:S) {
          z[i_10, i_11] ~ dbern(psi[i_10, i_11])
      }
  }
  for (i_12 in 1:M) {
      for (i_13 in 1:J) {
          for (i_14 in 1:S) {
              y[i_12, i_13, i_14] ~ dbern(p[i_12, i_13, i_14] *
                z[i_12, i_14])
          }
      }
  }
})

# Add missing values
y[1,1,] <- NA
y[2,2:5,] <- NA
y[3,,] <- NA

test_that("latent multispecies model works", {

  code <- nimbleCode({
    y[1:M,1:J,1:S] ~ MULTISPECIESOCCUPANCY(~x1[1:M], ~x2[1:M,1:J])
  })

  mod <- nimbleModel(code, constants = constants)

  expect_equal(mod$getCode(), code_ref)

  expect_equal(mod$z, apply(constants$y, c(1,3), max, na.rm=TRUE))
})

test_that("marginalized multispecies model works", {

  code <- nimbleCode({
    y[1:M,1:J,1:S] ~ MULTISPECIESOCCUPANCY(~x1[1:M], ~x2[1:M,1:J], marginalized = TRUE)
  })

  code_ref <- quote({
    for (i_3 in 1:M) {
        for (i_4 in 1:S) {
            logit(psi[i_3, i_4]) <- state_speciesID[speciesID[i_4]] +
                state_speciesID_x1[speciesID[i_4]] * x1[i_3]
        }
    }
    state_Intercept ~ dunif(-10, 10)
    state_x1 ~ dlogis(0, 1)
    state_sd_speciesID ~ dunif(0, 5)
    for (i_5 in 1:10) {
        state_speciesID[i_5] ~ dnorm(state_Intercept, sd = state_sd_speciesID)
    }
    state_sd_speciesID_x1 ~ dunif(0, 5)
    for (i_6 in 1:10) {
        state_speciesID_x1[i_6] ~ dnorm(state_x1, sd = state_sd_speciesID_x1)
    }
    for (i_7 in 1:M) {
        for (i_8 in 1:J) {
            for (i_9 in 1:S) {
                logit(p[i_7, i_8, i_9]) <- det_speciesID[speciesID[i_9]] +
                  det_speciesID_x2[speciesID[i_9]] * x2[i_7,
                    i_8]
            }
        }
    }
    det_Intercept ~ dunif(-10, 10)
    det_x2 ~ dlogis(0, 1)
    det_sd_speciesID ~ dunif(0, 5)
    for (i_10 in 1:10) {
        det_speciesID[i_10] ~ dnorm(det_Intercept, sd = det_sd_speciesID)
    }
    det_sd_speciesID_x2 ~ dunif(0, 5)
    for (i_11 in 1:10) {
        det_speciesID_x2[i_11] ~ dnorm(det_x2, sd = det_sd_speciesID_x2)
    }
    for (i_1 in 1:M) {
        for (i_2 in 1:S) {
            y[i_1, 1:J, i_2] ~ dOcc_v(probOcc = psi[i_1, i_2],
                probDetect = p[i_1, 1:J, i_2], len = J)
        }
    }
  })

  mod <- nimbleModel(code, constants = constants)

  expect_equal(mod$getCode(), code_ref)

})

test_that("nimbleOccu works with multispecies data",{
  sc <- list(x1 = x1)
  oc <- list(x2 = x2)
  mod <- nimbleOccu(~x1, ~x2, y, siteCovs=sc, obsCovs=oc, returnModel=TRUE)

  expect_equal(mod$getCode(), code_ref)
  
  # Species-level covariate
  sp <- list(x3 = rnorm(S), x4 = rnorm(3))

  mod <- nimbleOccu(~x1 + x3, ~x2, y, siteCovs=sc, obsCovs=oc,
                    speciesCovs = sp, returnModel=TRUE)
  
  expect_equal(mod$getCode()[[2]],
    quote({
      for (i_1 in 1:M) {
        for (i_2 in 1:S) {
        logit(psi[i_1, i_2]) <- state_x3 * x3[i_2] + state_speciesID[speciesID[i_2]] +
            state_speciesID_x1[speciesID[i_2]] * x1[i_1]
      }
    } 
    })[[2]]
  )

  mod <- nimbleOccu(~x1, ~x2 + x3, y, siteCovs=sc, obsCovs=oc,
                    speciesCovs = sp, returnModel=TRUE)
  
  expect_error(nimbleOccu(~x1, ~x2 + x4, y, siteCovs=sc, obsCovs=oc,
                    speciesCovs = sp, returnModel=TRUE))
})

# after, not sure if this is needed
nimbleOptions(enableModelMacros = macro_setting)
nimbleOptions(enableMacroComments = comm_setting)
