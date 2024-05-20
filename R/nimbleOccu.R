#' @export
nimbleOccu <- function(stateformula, detformula,
                       y, siteCovs = NULL, obsCovs = NULL,
                       statePriors=setPriors(intercept="dunif(-10, 10)", coefficient="dlogis(0, 1)", sd = "dunif(0,5)"), 
                       detPriors=setPriors(intercept="dunif(-10, 10)", coefficient="dlogis(0, 1)", sd = "dunif(0, 5)"),
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

  stopifnot(is.array(y))
  stopifnot(length(dim(y)) == 2 | length(dim(y)) == 3)
  stopifnot(ncol(y) > 1)

  M <- nrow(y)
  J <- ncol(y)
  
  S <- 1 # number of species
  if(length(dim(y)) == 3){
    S <- dim(y)[3]
    if(S == 1){
      # Drop last dimension if only one species
      y <- y[,,1,drop=FALSE]
    }
  }

  # Possibly step here to move data to front of vectors

  # Make data list
  data <- list()

  # Add z if latent model
  if(!marginalized){
    if(S == 1){
      z <- apply(y, 1, max, na.rm=TRUE)
    } else {
      z <- apply(y, c(1,3), max, na.rm=TRUE)
    }
    z[z==0] <- NA
    data$z <- z
  }

  # Make basic constants list
  constants <- list(y = y, M = M, J = J)  
  if(S > 1) constants <- c(constants, list(S = S))
  
  # Process site model and add required constants
  state_vars <- all.vars(stateformula)
  stopifnot(all(state_vars %in% names(siteCovs)))
  for (i in state_vars){
    stateformula <- addBracketToFormula(stateformula, str2lang(i), "[1:M]")
  }
  site_sub <- siteCovs[state_vars]
  # TODO: Check factors here?
  check_site <- sapply(site_sub, function(x){
    is.vector(x) & length(x) == M
  })
  stopifnot(all(check_site))
  constants <- c(constants, site_sub)

  # Process obs covariates
  obs_vars <- all.vars(detformula)
  stopifnot(all(obs_vars %in% c(names(obsCovs), names(siteCovs))))
  
  # Find covariates in the detection model that are actually site-level
  # modify formula and add to constants
  obs_in_site <- obs_vars[obs_vars %in% names(siteCovs)]
  obs_in_site <- obs_in_site[!obs_in_site %in% state_vars]
  if(length(obs_in_site) > 0){
    obs_site_sub <- siteCovs[obs_in_site]
    # Add required bracket ([1:M]) to the formula components
    for (i in obs_in_site){
      detformula <- addBracketToFormula(detformula, str2lang(i), "[1:M]")
    }
    check_site_obs <- sapply(obs_site_sub, function(x){
      is.vector(x) & length(x) == M
    })
    stopifnot(all(check_site_obs))
    constants <- c(constants, obs_site_sub)
  }

  # Find covariates in the detection model that are observation-level
  # modify formula and add to constants
  obs_vars <- obs_vars[!obs_vars %in% names(siteCovs)]
  if(length(obs_vars) > 0){
    obs_sub <- obsCovs[obs_vars]
    # Add required bracket ([1:M, 1:J])
    # TODO: Handle variable sample size for marginalized model (handling missing values)
    for (i in obs_vars){
      detformula <- addBracketToFormula(detformula, str2lang(i), "[1:M, 1:J]")
    }
    check_obs <- sapply(obs_sub[names(obs_sub) %in% names(obsCovs)], function(x){
      is.matrix(x) & all(dim(x) == c(M,J))
    })
    stopifnot(all(check_obs))
    constants <- c(constants, obs_sub)
  }
 
  # Generate code
  if(S == 1){
    # single species case
    resp <- quote(y[1:M, 1:J])
    macro <- quote(occupancy)
  } else {
    # multispecies
    resp <- quote(y[1:M, 1:J, 1:S])
    macro <- quote(multispeciesOccupancy)
  }

  code <- substitute({
    RESP ~ MACRO(stateformula = SF, detformula = DF,
                 statePriors = SP, detPriors = DP,
                 marginalized = MAR)
  }, list(RESP = resp, MACRO = macro, SF=stateformula, DF=detformula, 
          SP = statePriors, DP = detPriors, MAR = marginalized))

  args <- list(...)
  if(is.null(args$inits)) args$inits <- list()
  mod <- nimbleModel(code = code, data = data, constants = constants, 
                     inits = args$inits)
  
  if(returnModel){
    return(mod)
  }

  # Set up sampler
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

# Add dimension bracket (bracket; e.g. [1:M]) to a specific variable (target)
# found in a formula (inp)
addBracketToFormula <- function(inp, target, bracket){
  out <- addBracketToFormulaRecursive(inp, target, bracket)
  as.formula(as.call(out))
}

addBracketToFormulaRecursive <- function(inp, target, bracket){
  if(is.name(inp)){    
    if(inp == target){
      out <- str2lang(paste0(deparse(inp), bracket))
    } else {
      out <- inp
    }
  } else if(is.call(inp)) {
    out <- lapply(inp, addBracketToFormulaRecursive, target = target, bracket = bracket)
  } else {
    out <- inp
  }
  if(is.list(out)) out <- as.call(out)
  out
}
