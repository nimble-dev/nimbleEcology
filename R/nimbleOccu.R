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
#' @param samplerControl An optional list of control arguments to sampler functions 
#' to be passed `configureMCMC`. In addition, one can provide the element
#' \code{noncentered = TRUE} to add NIMBLE's "noncentered" sampler to jointly sample
#' each hyperparameter with its dependent random effects.
#' @param buildDerivs Build derivative capabilities for the model?
#' @param samplesOnly Save only the matrix of samples? If \code{FALSE}, the model is also saved.
#' @param ... Arguments passed to \code{runMCMC} such as \code{nchains}.
#'
#' @return Either MCMC samples, or a nimble model object if \code{returnModel = TRUE}.
#'
#' @export
nimbleOccu <- function(stateformula, detformula,
                       y, siteCovs = NULL, obsCovs = NULL, speciesCovs = NULL,
                       statePriors=setPriors(intercept="dunif(-10, 10)", coefficient="dlogis(0, 1)", sd = "dhalfflat()"), 
                       detPriors=setPriors(intercept="dunif(-10, 10)", coefficient="dlogis(0, 1)", sd = "dhalfflat()"),
                       marginalized = FALSE,
                       returnModel = FALSE,
                       sampler = c("default", "polyagamma", "RW_block", "barker", "hmc"),
                       savePsi = TRUE, saveP = FALSE,
                       samplerControl = list(),
                       buildDerivs = FALSE,
                       samplesOnly = FALSE,
                       ...){

  check_nimbleMacros_installed()

  ## TODO: decide on the actual default for the sampling scheme. 
  sampler <- match.arg(sampler)
  if(!returnModel) {  
    if(sampler %in% c("hmc", "barker"))
        buildDerivs <- TRUE
    ## Using block_RW and barker with latent model is conceptually fine but hasn't mixed well
    ## in CDFW examples.
    if(!marginalized && sampler %in% c("hmc", "block_RW", "barker")) {
      marginalized <- TRUE
      if(nimbleOptions("verbose")){
        message("  [Note] Using marginalized model as the ", sampler, " sampling is set up for use with the marginalized model.")
      }
    }
    if(marginalized && sampler == "polyagamma") {
      marginalized <- FALSE
      if(nimbleOptions("verbose")){
        message("  [Note] Using non-marginalized model as the ", sampler, " sampling is set up for use with the non-marginalized model.")
      }
    }
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
  if(!all(check_site)){
    wrong_dims <- names(check_site[!check_site])
    stop("Covariate(s) ", paste(wrong_dims, collapse=", "),
          " have incorrect dimensions", call.=FALSE)
  }
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
  if(!all(check_sp)){
    wrong_dims <- names(check_sp[!check_sp])
    stop("Covariate(s) ", paste(wrong_dims, collapse=", "),
          " have incorrect dimensions", call.=FALSE)
  }
  constants <- c(constants, sp_sub)

  # Process obs covariates
  obs_vars <- all.vars(detformula)
  stopifnot(all(obs_vars %in% c(names(obsCovs), names(siteCovs), names(speciesCovs))))
  
  # Find covariates in the detection model that are actually site-level
  # modify formula and add to constants
  obs_in_site <- obs_vars[obs_vars %in% names(siteCovs)]
  if(length(obs_in_site) > 0){
    obs_site_sub <- siteCovs[obs_in_site]
    # Add required bracket ([1:M]) to the formula components
    for (i in obs_in_site){
      detformula <- addBracketToFormula(detformula, str2lang(i), "[1:M]")
    }
    check_site_obs <- sapply(obs_site_sub, function(x){
      is.vector(x) & length(x) == M
    })
    if(!all(check_site_obs)){
      wrong_dims <- names(check_site_obs[!check_site_obs])
      stop("Covariate(s) ", paste(wrong_dims, collapse=", "),
           " have incorrect dimensions", call.=FALSE)
    }
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
    if(!all(check_sp_obs)){
      wrong_dims <- names(check_sp_obs[!check_sp_obs])
      stop("Covariate(s) ", paste(wrong_dims, collapse=", "),
           " have incorrect dimensions", call.=FALSE)
    }
    obs_sp_sub <- obs_sp_sub[!names(obs_sp_sub) %in% names(constants)]
    constants <- c(constants, obs_sp_sub)
  }

  # Find covariates in the detection model that are observation-level
  # modify formula and add to constants
  obs_vars <- obs_vars[obs_vars %in% names(obsCovs)]
  if(length(obs_vars) > 0){
    obs_sub <- obsCovs[obs_vars]
    # Add required bracket ([1:M, 1:J])
    for (i in obs_vars){
      detformula <- addBracketToFormula(detformula, str2lang(i), "[1:M, 1:J]")
    }
    check_obs <- sapply(obs_sub[names(obs_sub) %in% names(obsCovs)], function(x){
      is.matrix(x) & all(dim(x) == c(M,J))
    })
    if(!all(check_obs)){
      wrong_dims <- names(check_obs[!check_obs])
      stop("Covariate(s) ", paste(wrong_dims, collapse=", "),
           " have incorrect dimensions", call.=FALSE)
    }
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
                     inits = args$inits, buildDerivs = buildDerivs)
  
  if(returnModel){
    return(mod)
  }

  conf <- configureOccuMCMC(mod, sampler, samplerControl, S = S,
                            marginalized = marginalized, occ_var = 'z')
 
  if(savePsi){
    conf$addMonitors("psi")
  }
  if(saveP){
    conf$addMonitors("p")
  }

  mcmc <- buildMCMC(conf, print = nimbleOptions("verbose"))
  
  modC <- compileNimble(mod)
  mcmcC <- compileNimble(mcmc, project = mod)

  samples <- do.call(runMCMC, c(list(mcmc = mcmcC), args))
  if(samplesOnly) return(samples)

  out <- list(model = mod, samples = samples, call = match.call())
  class(out) <- c("nimbleOccu", class(out))
  out
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

configureOccuMCMC <- function(model, sampler, samplerControl, S, marginalized, occ_var = 'z') {

    ## TODO: consider site-level covariates and obs-level:
    ## how do we set up blocking when these are present.
    ## Possibly we would add each of these to each species-level block, since
    ## they affect likelihood for all species. 
    
    ## TODO: More generally, this code basically assumes model structure like the CDFW examples.
    ## We need to check it for other possible specifications, including cases
    ## where effects do not vary by species.

    # Get all parameters created by macros
    macro_pars <- model$getMacroParameters()
  
    # Extract just coefficients (which is all LINPRED will create)
    coefs <- macro_pars[["nimbleMacros::LINPRED"]]
    stopifnot(length(coefs) == 2) # possibly unnecessary checks
    state_coeffs_vars <- unique(coefs[[1]])
    stopifnot(all(grepl("^state_", state_coeffs_vars)))
    detect_coeffs_vars <- unique(coefs[[2]])
    stopifnot(all(grepl("^det_", detect_coeffs_vars)))

    # Get hyperparameters
    # Start by getting all coefficients and hyperparameters from the priors macros
    # No need to distinguish between state/det here.
    coefs_hyperpars <- unlist(macro_pars[["nimbleMacros::LINPRED_PRIORS"]])
    # Remove the coefficients to leave only hyperparameters
    hyper_vars <- coefs_hyperpars[! coefs_hyperpars %in% unlist(coefs)]
    hyper_vars <- unique(hyper_vars) # just in case

    # This is the only way I can think of to split location and scale.
    hyper_vars_scale <- hyper_vars[grep("_sd_", hyper_vars)]
    hyper_vars_location <- hyper_vars[grep("_sd_", hyper_vars, invert = TRUE)]
    if(length(setdiff(hyper_vars, c(hyper_vars_scale, hyper_vars_location))))
        stop("Detected unexpected hyperparameters that are neither location nor scale parameters.")

    if(S > 1){
        # Extract only coefs random by species (needed for block samplers)
        # all coefs random by species will contain _speciesID
        state_vars_sp <- state_coeffs_vars[grepl("_speciesID", state_coeffs_vars)]
        detect_vars_sp <- detect_coeffs_vars[grepl("_speciesID", detect_coeffs_vars)]
        vars_sp <- c(state_vars_sp, detect_vars_sp)

        # List of coefs that are not random by species
        state_vars_notsp <- state_coeffs_vars[!grepl("_speciesID", state_coeffs_vars)]
        detect_vars_notsp <- detect_coeffs_vars[!grepl("_speciesID", detect_coeffs_vars)]
        vars_notsp <- c(state_vars_notsp, detect_vars_notsp)
    }

    if(sampler == "hmc") {
        conf <- configureMCMC(model, nodes = NULL, print = FALSE)
        conf$addSampler(nimbleHMC::sampler_NUTS, target = conf$getUnsampledNodes(),
          control = list(warmupMode = "default"))
    }

    # Note that in the species-specific block sampling below,
    # we include non-species-specific covariate coeffs within each block,
    # since they are part of the linear predictor for each species.
    # They will be sampled repeatedly.

    
    if(sampler == "polyagamma") {
        # Configure hyperparameters if they exist.
        conf <- configureMCMC(model, nodes = hyper_vars, print = FALSE)
        # Add configuration for latent state.
        conf$addSampler(occ_var, "jointBinary", control = list(order = TRUE))

        if(S == 1){
            # For single species, just one sampler each for state/detection
            conf$addSampler(type = "sampler_polyagamma", target= state_coeffs_vars, 
                            control = list(fixedDesignColumns = TRUE))
            conf$addSampler(type = "sampler_polyagamma", target= detect_coeffs_vars, 
                            control = list(fixedDesignColumns = TRUE))
        } else {
            for (i in 1:S){
                pgNodes <- paste0(state_vars_sp, "[", i, "]")
                pgNodes <- c(pgNodes, state_vars_notsp) 
                conf$addSampler(type = "sampler_polyagamma", target=pgNodes, 
                                control = list(fixedDesignColumns = TRUE))
                pgNodes <- paste0(detect_vars_sp, "[", i, "]")
                pgNodes <- c(pgNodes, detect_vars_notsp)
                conf$addSampler(type = "sampler_polyagamma", target=pgNodes, 
                                control = list(fixedDesignColumns = TRUE))
            }
        }
    }
    
    if(sampler == "default") {
        conf <- configureMCMC(model, 
                              nodes = c(hyper_vars, detect_coeffs_vars, state_coeffs_vars), 
                              print = FALSE)
        if(!marginalized)
            conf$addSampler(occ_var, "jointBinary", control = list(order = TRUE))
    }
    
    if(sampler == "RW_block") {
        # Non-species-specific parameters included here
        conf <- configureMCMC(model, nodes = hyper_vars, print = FALSE)
        if(S == 1) {
            conf$addSampler(c(state_coeffs_vars, detect_coeffs_vars), type = "RW_block",
                            print=FALSE, control = samplerControl)

        } else {
            ## Species-specific slopes and intercepts
            for (i in 1:S){
                blockNodes <- c(vars_notsp, paste0(vars_sp, "[", i, "]"))
                conf$addSampler(blockNodes, type = "RW_block", print=FALSE,
                                control = samplerControl)
            }
        }
        if(!marginalized)  # At the moment, this is ruled out in initial `nimbleOccu` checks.
          conf$addSampler(occ_var, "jointBinary", control = list(order = TRUE))
    }
    
    if(sampler == "barker") {
        conf <- configureMCMC(model, nodes = hyper_vars, print = FALSE)
        if(S == 1){
            conf$addSampler(state_coeffs_vars, type = "barker", print=FALSE, control = samplerControl)
            conf$addSampler(detect_coeffs_vars, type = "barker", print=FALSE, control = samplerControl)
        } else {
            ## Species-specific slopes and intercepts
            for (i in 1:S){
              blockNodes <- c(state_vars_notsp, paste0(state_vars_sp, "[", i, "]"))
              conf$addSampler(blockNodes, type = "barker", print=FALSE, control = samplerControl)
              blockNodes <- c(detect_vars_notsp, paste0(detect_vars_sp, "[", i, "]"))
              conf$addSampler(blockNodes, type = "barker", print=FALSE, control = samplerControl)
            }
        }
        if(!marginalized)  # At the moment, this is ruled out in initial `nimbleOccu` checks.
            conf$addSampler(occ_var, "jointBinary", control = list(order = TRUE))
    }

    if(sampler != "hmc"){
        ## Use slice sampling for standard deviations if not conjugate.
        for(node in model$expandNodeNames(hyper_vars_scale)) {
            samplers <- conf$getSamplers(node)
            if(length(samplers) == 1 && length(grep("^conjugate", samplers[[1]]$name, invert = TRUE))) 
                conf$replaceSamplers(node, "slice")
        }
        
        
        ## Optionally add noncentered sampling for each hyperparameter with its random effects.
        noncentered <- samplerControl$noncentered
        if(is.null(noncentered)) noncentered <- FALSE
        if(noncentered) {
            for(param in conf$model$expandNodeNames(hyper_vars_location))
                conf$addSampler(type = 'noncentered', target = param,
                                control = list(sampler = "RW", param = "location", log = TRUE))
            for(param in conf$model$expandNodeNames(hyper_vars_scale))
                conf$addSampler(type = 'noncentered', target = param,
                                control = list(sampler = "RW", param = "scale", log = TRUE))
        }
    }
    return(conf)
}





## "Vectorized" occupancy distribution over multiple sites.
## `dOcc_v` in `nimbleEcology only vectorizes over visits at a single site.
dOcc_v_multisite <- nimbleFunction(run = function(x = double(2),
                                            probOcc = double(1), probDetect = double(2),
                                            lens = double(1), n.sites = double(0),
                                            log = integer(0, default = 0)) {
    returnType(double(0))
    
    n.sites_noDeriv <- 0L 
    len <- 0L 
    n.sites_noDeriv <- ADbreak(n.sites)

    logProb <- 0
    for(j in 1:n.sites_noDeriv) {
        len  <- ADbreak(lens[j])
        logProb <- logProb + dOcc_v(x[j, 1:len], probOcc = probOcc[j], probDetect = probDetect[j, 1:len], len = len, log=TRUE)
    }
    if(log) return(logProb) else return(exp(logProb))
}, buildDerivs = list(run = list(ignore=c('j')))
)

## This function causes 'package 'methods' is used but not declared'
## is it just because it literally contains the word "methods"?
## "Vectorized" sampler that samples all occupancy indicators at once.
#' @export
sampler_jointBinary <- nimbleFunction(
    name = 'sampler_jointBinary',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        targetAsScalar <- targetAsScalar[!model$isData(targetAsScalar)]

        
        ## Put in same order as configureMCMC would do so that sampling exactly matches
        ## order of application of binary samplers.
        if(!is.null(control$order) && control$order) {
            fullNodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
            mtch <- match(targetAsScalar, fullNodes)
            targetAsScalar <- targetAsScalar[order(mtch)]
        }            
            
        ccList <- nimble:::mcmc_determineCalcAndCopyNodes(model, targetAsScalar)
        calcNodes <- ccList$calcNodes; calcNodesNoSelf <- ccList$calcNodesNoSelf; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch
        p <- length(targetAsScalar)

        lens <- rep(0, p)
        first <- 1
        for(i in 1:p) {
            calcNodesOne <- model$getDependencies(targetAsScalar[i], self = TRUE)  
            lens[i] <- length(calcNodesOne)
            last <- first+lens[i]-1
            calcNodes[first:last] <- calcNodesOne
            first <- last + 1
        }
        currentLogProb <- otherLogProb <- rep(0, p)
    },
    run = function() {
        start <- 1
        for(i in 1:p) {
            end <- start + lens[i] - 1
            currentLogProb[i] <<- model$getLogProb(calcNodes[start:end])
            start <- end + 1
        }
        values(model, targetAsScalar) <<- 1 - values(model, targetAsScalar)
        model$calculate(calcNodes)
        start <- 1
        for(i in 1:p) {
            end <- start + lens[i] - 1
            otherLogProb[i] <<- model$getLogProb(calcNodes[start:end])
            start <- end + 1
        }
        acceptanceProb <- 1/(exp(currentLogProb - otherLogProb) + 1)
        jump <- (!is.nan(acceptanceProb)) & (runif(p,0,1) < acceptanceProb)

        for(i in 1:p) 
            if(!jump[i]) 
                values(model, targetAsScalar[i]) <<- 1-values(model, targetAsScalar[i])

        model$calculate(calcNodes)
        
        nimCopy(from = model, to = mvSaved, row = 1, nodes = targetAsScalar, logProb = TRUE)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
    },
    methods = list(
        reset = function() { }
    )
)

# nimbleOccu methods

#' @export
summary.nimbleOccu <- function(object, levels=c(0.025, 0.5, 0.975), digits=2, ...){

  # Convert to coda object
  if(is.list(object$samples)){
    object$samples <- coda::as.mcmc.list(lapply(object$samples, coda::as.mcmc))
    pars <- colnames(object$samples[[1]])
  } else {
    object$samples <- coda::as.mcmc(object$samples)
    pars <- colnames(object$samples)
  }

  state_tab <- summary_table(pars[grepl("^state_", pars)], object$samples, levels)
  det_tab <- summary_table(pars[grepl("^det_", pars)], object$samples, levels)

  cat("Occupancy:\n")
  print(round(state_tab, digits))
  cat("\nDetection:\n")
  print(round(det_tab, digits))

  invisible(list(state=state_tab, det=det_tab))
}

summary_table <- function(pars, samples, levels){
  sumtab <- lapply(pars, function(x){
    sampx <- unlist(samples[,x])
    out <- c(mean = mean(sampx, na.rm=TRUE),
    stats::quantile(sampx, probs=levels, na.rm=TRUE))
    if(is.list(samples)){
      rhat <- coda::gelman.diag(samples[,x])
      out <- c(out, Rhat=rhat$psrf[1])
    }
    out <- c(out, effSize=unname(coda::effectiveSize(samples[,x])))
    out
  })
  sumtab <- as.data.frame(do.call(rbind, sumtab))
  rownames(sumtab) <- gsub("^state_|^det_", "", pars)
  sumtab
}

#' @export
print.nimbleOccu <- function(x, ...){
  summary(x, ...)
}

