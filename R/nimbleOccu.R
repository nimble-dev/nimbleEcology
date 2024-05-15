#' @export
nimbleOccu <- function(stateformula, detformula,
                       y, siteCovs = NULL, obsCovs = NULL, 
                       marginalized = FALSE,
                       returnModel = FALSE,
                       sampler = c("default", "hmc"),
                       savePsi = TRUE, saveP = FALSE,
                       ...){

  check_nimbleMacros_installed()
  
  sampler <- match.arg(sampler)
  if(sampler == "hmc"){
    stop("Doesn't work until nimbleEcology AD branch is merged", call.=FALSE)
    if(!marginalized){
      stop("hmc sampler requires marginalized model", call.=FALSE)
    }
    nimbleOptions(buildModelDerivs=TRUE) # this is a bit hacky
  }

  stopifnot(is.matrix(y))
  stopifnot(ncol(y) > 1)

  M <- nrow(y)
  J <- ncol(y)

  # Possibly step here to move data to front of vectors

  # Make data list
  data <- list()
  if(!marginalized){
    z <- apply(y, 1, max, na.rm=TRUE)
    z[z==0] <- NA
    data$z <- z
  }

  # Make constants list
  constants <- list(y = y, M = M, J = J)
  
  # Site covariates
  state_vars <- all.vars(stateformula)
  stopifnot(all(state_vars %in% names(siteCovs)))
  site_sub <- siteCovs[state_vars]
  # Check factors here?
  check_site <- sapply(site_sub, function(x){
    is.vector(x) & length(x) == M
  })
  stopifnot(all(check_site))
  constants <- c(constants, site_sub)

  # Obs covariates
  obs_vars <- all.vars(detformula)
  stopifnot(all(obs_vars %in% names(obsCovs)))
  obs_sub <- obsCovs[obs_vars]
  check_obs <- sapply(obs_sub, function(x){
    is.matrix(x) & all(dim(x) == c(M,J))
  })
  stopifnot(all(check_obs))
  constants <- c(constants, obs_sub)
  
  code <- substitute({
    y[1:M, 1:J] ~ occupancy(stateformula = SF,
                            detformula = DF,
                            marginalized = MAR)
  }, list(SF=stateformula, DF=detformula, MAR = marginalized))

  args <- list(...)
  mod <- nimbleModel(code = code, data = data, constants = constants, 
                     inits = args$inits)
  
  if(returnModel){
    return(mod)
  }

  if(sampler == "default"){
    conf <- configureMCMC(mod, print = FALSE)
  } else if(sampler == "hmc"){
    if(!requireNamespace("nimbleHMC")){
      stop("Install nimbleHMC package", call.=FALSE)
    }
    if(is.null(args$warmup)){
      stop("Must provide warmup argument", call.=FALSE)
    }
    conf <- nimbleHMC::configureHMC(mod, print = FALSE,
              control = list(warmupMode = "iterations", warmup=args$warmup))
  }

  if(savePsi){
    conf$addMonitors("psi")
  }
  if(saveP){
    conf$addMonitors("p")
  }

  mcmc <- buildMCMC(conf)
  
  modC <- compileNimble(mod)
  mcmcC <- compileNimble(mcmc, project = mod)

  runMCMC(mcmcC, ...)
}
