#' Fit occupancy models with nimble
#'
#' Fit single or multi-species occupancy models with nimble. The state
#' and detection models are defined with R formulas and data are supplied to
#' the y, siteCovs, obsCovs, and speciesCovs arguments.
#' You may either fit the model or return the model object for further processing..
#'
#' @name nimbleOccu
#' @author Ken Kellner
#'
#' @param stateformula An R formula for the occupancy/state model, possibly with the 
#'  parameters followed by brackets containing indices.
#' @param detformula An R formula for the detection model, possibly with the 
#'  parameters followed by brackets containing indices.
#' @param y The detection/non-detection data. For single-species models this 
#'  should be an MxJ matrix; for multi-species models it should be an MxJxS
#'  array, where M is sites, J is occasions, and S is species.
#' @param siteCovs A named list of site-level covariates. Each should be a vector of length M.
#' @param obsCovs A named list of observation-level covariates. Each should be 
#'  a matrix of dimensions MxJ.
#' @param speciesCovs A named list of species-level covariates. Each should be
#'  a vector of length S. Only used in multi-species model.
#' @param statePriors Prior specifications for state model parameters, should be 
#'  generated with nimbleMacros::setPriors()
#' @param detPriors Prior specifications for det model parameters, should be generated 
#'  with nimbleMacros::setPriors()
#' @param marginalized Logical. If TRUE, fit the marginalized model using the
#'  dOcc_v nimbleFunction
#' @param returnModel Logical. If TRUE, the nimble model object is returned.
#' @param sampler MCMC sampler to use. Currently only the nimble default sampler works.
#' @param savePsi Save site-level occupancy probabilities in the output?
#' @param saveP Save site x occasion level detection probabilities in the output?
#' @param ... Arguments passed to \code{runMCMC} such as \code{nchains}.
#'
#' @return Either MCMC samples, or a nimble model object if \code{returnModel = TRUE}.
#'
#' @export
nimbleOccu <- function(stateformula, detformula,
                       y, siteCovs = NULL, obsCovs = NULL, speciesCovs = NULL,
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
    if(!marginalized){
      message("setting marginalized = TRUE")
      marginalized <- TRUE
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

  if(S == 1 & !is.null(speciesCovs)){
    stop("Can't have species covariates in a single-species model", call.=FALSE)
  }

  # Possibly step here to move data to front of vectors

  # Make data list
  data <- list()

  # Add z if latent model
  if(!marginalized){
    if(S == 1){
      z <- apply(y, 1, function(x) if(all(is.na(x))) return(0) else max(x, na.rm=TRUE))
    } else {
      z <- apply(y, c(1,3), function(x) if(all(is.na(x))) return(0) else max(x, na.rm=TRUE))
    }
    z[z==0] <- NA
    data$z <- z
  }

  # Make basic constants list
  constants <- list(y = y, M = M, J = J)  
  if(S > 1) constants <- c(constants, list(S = S))
  
  # Process site model and add required constants
  state_vars <- all.vars(stateformula)
  stopifnot(all(state_vars %in% c(names(siteCovs), names(speciesCovs))))

  state_vars_sc <- state_vars[state_vars %in% names(siteCovs)]
  for (i in state_vars_sc){
    stateformula <- addBracketToFormula(stateformula, str2lang(i), "[1:M]")
  }
  site_sub <- siteCovs[state_vars_sc]
  # TODO: Check factors here?
  check_site <- sapply(site_sub, function(x){
    is.vector(x) & length(x) == M
  })
  stopifnot(all(check_site))
  constants <- c(constants, site_sub)

  state_vars_sp <- state_vars[state_vars %in% names(speciesCovs)]
  for (i in state_vars_sp){
    stateformula <- addBracketToFormula(stateformula, str2lang(i), "[1:S]")
  }
  sp_sub <- speciesCovs[state_vars_sp]
  # TODO: Check factors here?
  check_sp <- sapply(sp_sub, function(x){
    is.vector(x) & length(x) == S
  })
  stopifnot(all(check_sp))
  constants <- c(constants, sp_sub)

  # Process obs covariates
  obs_vars <- all.vars(detformula)
  stopifnot(all(obs_vars %in% c(names(obsCovs), names(siteCovs), names(speciesCovs))))
  
  # Find covariates in the detection model that are actually site-level
  # modify formula and add to constants
  obs_in_site <- obs_vars[obs_vars %in% names(siteCovs)]
  #obs_in_site <- obs_in_site[!obs_in_site %in% state_vars] # TODO: I think this check is wrong
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
    obs_site_sub <- obs_site_sub[!names(obs_site_sub) %in% names(constants)]
    constants <- c(constants, obs_site_sub)
  }

  # Find covariates in the detection model that are actually species-level
  # modify formula and add to constants
  obs_in_sp <- obs_vars[obs_vars %in% names(speciesCovs)]
  if(length(obs_in_sp) > 0){
    obs_sp_sub <- speciesCovs[obs_in_sp]
    # Add required bracket ([1:S]) to the formula components
    for (i in obs_in_sp){
      detformula <- addBracketToFormula(detformula, str2lang(i), "[1:S]")
    }
    check_sp_obs <- sapply(obs_sp_sub, function(x){
      is.vector(x) & length(x) == S
    })
    stopifnot(all(check_sp_obs))
    obs_sp_sub <- obs_sp_sub[!names(obs_sp_sub) %in% names(constants)]
    constants <- c(constants, obs_sp_sub)
  }

  # Find covariates in the detection model that are observation-level
  # modify formula and add to constants
  obs_vars <- obs_vars[obs_vars %in% names(obsCovs)]
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
    macro <- quote(OCCUPANCY)
    code <- substitute({
    RESP ~ OCCUPANCY(stateformula = SF, detformula = DF,
                 statePriors = SP, detPriors = DP,
                 marginalized = MAR)
    }, list(RESP = resp, SF=stateformula, DF=detformula, 
          SP = statePriors, DP = detPriors, MAR = marginalized))
  } else {
    # multispecies
    resp <- quote(y[1:M, 1:J, 1:S])
    spcovs <- ""
    if(!is.null(speciesCovs)) spcovs <- names(speciesCovs)
    code <- substitute({
    RESP ~ MULTISPECIESOCCUPANCY(stateformula = SF, detformula = DF,
                 statePriors = SP, detPriors = DP,
                 marginalized = MAR, speciesCovs=SPCOVS)
    }, list(RESP = resp, SF=stateformula, DF=detformula, 
          SP = statePriors, DP = detPriors, MAR = marginalized, SPCOVS=spcovs))
  }

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
    args$warmup <- NULL
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

  do.call(runMCMC, c(list(mcmc = mcmcC), args))
}

# Add dimension bracket (bracket; e.g. [1:M]) to a specific variable (target)
# found in a formula (inp)
addBracketToFormula <- function(inp, target, bracket){
  out <- addBracketToFormulaRecursive(inp, target, bracket)
  stats::as.formula(as.call(out))
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
