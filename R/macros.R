#' Macro to fit single-season, single-species occupancy model
#'
#' Generates nimble code for a single-season, single-species occupancy model.
#' Covariates can be specified for both the detection and occupancy (state)
#' submodels using R formulas.
#'
#' @name occupancy
#' @author Ken Kellner
#' 
#' @param stateformula An R formula for the occupancy/state model, possibly with the 
#'  parameters followed by brackets containing indices.
#' @param detformula An R formula for the detection model, possibly with the 
#'  parameters followed by brackets containing indices.
#' @param statePrefix All state model coefficient names will begin with this prefix.
#'  default is state_ (so x becomes state_x, etc.)
#' @param detPrefix All detection model coefficient names will begin with this prefix.
#'  default is det_ (so x becomes det_x, etc.)
#' @param statePriors Prior specifications for state model parameters, should be 
#'  generated with nimbleMacros::setPriors()
#' @param detPriors Prior specifications for det model parameters, should be generated 
#'  with nimbleMacros::setPriors()
#' @param marginalized Logical. If TRUE, fit the marginalized model using the
#'  dOcc_v nimbleFunction
#' @param centerVar Grouping covariate to 'center' on in parameterization. By
#'  default all random effects have mean 0 as with lme4.
#'
#' @details This macro should be placed on the right-hand side of an assignment
#'  (it doesn't matter if you use <- or ~). The left-hand side must be a
#'  matrix of detection-nondetection data with dimensions M x J, where
#'  M is sites and J is occasions. This matrix may contain
#'  missing values (e.g. if you have different numbers of occasions by site). However,
#'  to use the marginalized form of the model, the missing values for site M
#'  must come at the end of the row of occasion data 1:J. For example
#'  if you have J = 5 but the first site was only sampled 3 times, the 
#'  corresponding row of the matrix should be [1, 1, 1, NA, NA] rather than
#'  for example [1, NA, 1, 1, NA].
#'
#' @examples
#' nimbleOptions(enableModelMacros = TRUE)
#' y <- matrix(rbinom(10, 1, 0.5), 5, 2)
#' x <- rnorm(5)
#' const <- list(y=y, x=x, M=nrow(y), J = ncol(y))
#'
#' code <- nimbleCode({
#'   y[1:M, 1:J] ~ occupancy(~x[1:M], ~1)
#' })
#'
#' mod <- nimbleModel(code, constants=const)
#' mod$getCode()
NULL

#' @export
occupancy <- nimble::model_macro_builder(
function(stoch, LHS, stateformula = ~1, detformula = ~1,
         statePrefix=quote(state_), detPrefix=quote(det_), 
         statePriors=setPriors(intercept="dunif(-10, 10)", coefficient="dlogis(0, 1)"), 
         detPriors=setPriors(intercept="dunif(-10, 10)", coefficient="dlogis(0, 1)"),
         marginalized=FALSE, centerVar=NULL,
         modelInfo, .env){
    
    check_nimbleMacros_installed()

    site_dim <- LHS[[3]] # e.g. 1:M
    obs_dim <- LHS[[4]]  # e.g. 1:J
    
    if(marginalized){
      yraw <- eval(LHS[[2]], envir=modelInfo$constants)
      yrng <- eval(site_dim, envir=modelInfo$constants)
      nsamp <- apply(yraw, 1, function(x) sum(!is.na(x)))
      if(any(nsamp == 1)){
        split_sites <- TRUE
        modelInfo$constants$sidx1 <- yrng[which(nsamp > 1)]
        modelInfo$constants$M1 <- length(modelInfo$constants$sidx1)
        modelInfo$constants$sidx2 <- yrng[which(nsamp == 1)]
        modelInfo$constants$M2 <- length(modelInfo$constants$sidx2)
        # After sorting NAs this won't be needed (?)
        modelInfo$constants$vidx <- apply(yraw[nsamp==1,,drop=FALSE], 1, function(x) which(!is.na(x))) 
        single_site_loop <- substitute({
          Y[sidx2[1:M2], vidx[1:M2]] ~ nimbleMacros::forLoop(dbern(psi[sidx2[1:M2]] * p[sidx2[1:M2], vidx[1:M2]]))
        }, list(Y=LHS[[2]], OBSDIM=obs_dim))
      }
    }
    
    # Code to calculate parameters psi and p
    param_code <- substitute({
      psi[SITEDIM] <- nimbleMacros::linPred(STATEFORM, link=logit, coefPrefix=STATEPREFIX, 
                                sdPrefix=STATEPREFIX, priorSettings=STATEPRIOR, 
                                centerVar=CENTVAR)  
      p[SITEDIM, OBSDIM] <- nimbleMacros::linPred(DETFORM, link=logit, coefPrefix=DETPREFIX, 
                                      sdPrefix=DETPREFIX, priorSettings=DETPRIOR, 
                                      centerVar=CENTVAR)
      }, 
      list(SITEDIM=site_dim, OBSDIM=obs_dim,
           STATEFORM=stateformula, DETFORM = detformula,
           STATEPREFIX=statePrefix, DETPREFIX = detPrefix,
           STATEPRIOR=statePriors, DETPRIOR=detPriors, CENTVAR=centerVar)
    )
   
    # Code for y model
    if(!marginalized){
      # Latent model
      ymod_code <- substitute({  
        z[SITEDIM] ~ nimbleMacros::forLoop(dbern(psi[SITEDIM]))
        Y ~ nimbleMacros::forLoop(dbern(p[SITEDIM, OBSDIM]*z[SITEDIM]))
        }, 
        list(Y=LHS, SITEDIM=site_dim, OBSDIM=obs_dim)
      )

      # Add latent z to inits (ideally would be data but can't at the moment)
      yname <- deparse(LHS[[2]])
      if(!yname %in% names(modelInfo$constants)){
        stop(paste0("The detection-nondetection matrix ", yname, " should be in the constants"), 
            call.=FALSE)
      }
      y <- modelInfo$constants[[yname]]
      z <- apply(y, 1, max, na.rm=TRUE)
      #z[z == 0] <- NA
      #modelInfo$constants$z <- z
      modelInfo$inits <- list(z = z)

    } else {
      # Marginalized model
      # Because of the way the dOcc_v function works we can't use 
      # the forLoop macro, so have to write the loop code manually

      # New index variable for site loop
      site_idx <- as.name(modelInfo$indexCreator())
      site_idx2 <- site_idx

      if(split_sites){
        site_dim <- quote(1:M1)
        site_idx2 <- substitute(sidx1[IDX], list(IDX=site_idx))
      }

      # Get observation range (e.g. 1:J)
      obs_dim2 <- replaceIndex(obs_dim, site_dim, site_idx) 
      # Number of samples
      nsamples <- obs_dim2[[3]]

      ymod_code <- substitute({  
        for (SITEIDX in SITEDIM){
          Y[SITEIDX2, OBSDIM2] ~ dOcc_v(
            probOcc = psi[SITEIDX2],
            probDetect = p[SITEIDX2, OBSDIM2],
            len = NSAMP
            )
        }
        }, 
        list(Y=LHS[[2]], SITEDIM=site_dim, OBSDIM=obs_dim,
             SITEIDX=site_idx, SITEIDX2=site_idx2, OBSDIM2=obs_dim2, NSAMP=nsamples)
      )

      if(split_sites){
        ymod_code <- embedLinesInCurlyBrackets(list(ymod_code, single_site_loop))
      }
    }

    code <- embedLinesInCurlyBrackets(list(param_code, ymod_code))

    # Return code and model info
    list(code=code, modelInfo=modelInfo)
  },
use3pieces=TRUE,
unpackArgs=TRUE
)
class(occupancy) <- "model_macro"

setPriors <- function(...){
  check_nimbleMacros_installed()
  nimbleMacros::setPriors(...)
}

check_nimbleMacros_installed <- function(){
  check <- requireNamespace("nimbleMacros")
  if(!check){
    stop("Install nimbleMacros package", call.=FALSE)
  }
  #assign("forLoop", nimbleMacros::forLoop, envir = globalenv())
  #assign("linPred", nimbleMacros::linPred, envir = globalenv())
  #assign("priors", nimbleMacros::priors, envir = globalenv())
  invisible()
}

# From nimbleMacros
hasBracket <- function(code, recursive=TRUE){
  if (length(code) < 2) return(FALSE)
  if (code[[1]] == "[") return(TRUE)
  if(recursive){
    if (is.name(code[[1]])) return(hasBracket(code[[2]]))
  }
    FALSE
}

hasAdjustment <- function(code){
  if(is.name(code)) return(FALSE)
  if(length(code) < 3) return(FALSE)
  code[[1]] == "-" | code[[1]] == "+"
}

# Replace a provided index in some code with a new value
replaceIndex <- function(code, old_idx, new_idx){
  #stopifnot(hasBracket(code))
  code_list <- as.list(code)
  code_list <- lapply(code_list, function(x){
    if(hasBracket(x)){
      return(replaceIndex(x, old_idx, new_idx))
    }
    if(hasAdjustment(x)){
      return(replaceIndex(x, old_idx, new_idx))
    }
    x
  })
  idx <- which(code_list == old_idx)
  # If old index is not found do nothing
  if(length(idx) != 1) return(as.call(code_list))
  code_list[[idx]] <- new_idx
  as.call(code_list)
}

# Function copied from the nimble package
embedLinesInCurlyBrackets <- function(lines) {
  as.call(c(list(quote(`{`)), lines))
}


#' Macro to fit single-season, multi-species occupancy model
#'
#' Generates nimble code for a single-season, multi-species occupancy model,
#' in which all intercepts and slopes are assumed to be random effects by species.
#' Covariates can be specified for both the detection and occupancy (state)
#' submodels using R formulas.
#'
#' @name multispeciesOccupancy
#' @author Ken Kellner
#' 
#' @param stateformula An R formula for the occupancy/state model, possibly with the 
#'  parameters followed by brackets containing indices.
#' @param detformula An R formula for the detection model, possibly with the 
#'  parameters followed by brackets containing indices.
#' @param statePrefix All state model coefficient names will begin with this prefix.
#'  default is state_ (so x becomes state_x, etc.)
#' @param detPrefix All detection model coefficient names will begin with this prefix.
#'  default is det_ (so x becomes det_x, etc.)
#' @param statePriors Prior specifications for state model parameters, should be 
#'  generated with nimbleMacros::setPriors()
#' @param detPriors Prior specifications for det model parameters, should be generated 
#'  with nimbleMacros::setPriors()
#' @param marginalized Logical. If TRUE, fit the marginalized model using the
#'  dOcc_v nimbleFunction
#' @param speciesID Name to use for species ID index vector 
#'
#' @details This macro should be placed on the right-hand side of an assignment
#'  (it doesn't matter if you use <- or ~). The left-hand side must be an
#'  array of detection-nondetection data with dimensions M x J x S, where
#'  M is sites, J is occasions, and S is species.  This array may contain
#'  missing values (e.g. if you have different numbers of occasions by site). However,
#'  to use the marginalized form of the model, the missing values for site M
#'  and species S must come at the end of the occasion data 1:J. For example
#'  if you have J = 5 but the first site was only sampled 3 times, the 
#'  corresponding row of the array should be [1, 1, 1, NA, NA] rather than
#'  for example [1, NA, 1, 1, NA].
#'
#' @examples
#' nimbleOptions(enableModelMacros = TRUE)
#' M <- 100 # sites
#' J <- 5   # occasions
#' S <- 10  # species
#' y <- array(rbinom(M*J*S, 1, 0.5), c(M,J,S))
#' x <- rnorm(M) # site covariate
#'
#' const <- list(M=M, J=J, S=S, y=y, x=x)
#'
#' code <- nimbleCode({
#'   y[1:M, 1:J, 1:S] ~ multispeciesOccupancy(~x[1:M], ~1)
#' })
#'
#' mod <- nimbleModel(code, constants=const)
#' mod$getCode()
NULL

#' @export
multispeciesOccupancy <- nimble::model_macro_builder(
function(stoch, LHS, stateformula, detformula, 
         statePrefix=quote(state_), detPrefix=quote(det_), 
         statePriors=setPriors(intercept="dunif(-10, 10)", coefficient="dlogis(0, 1)", sd = "dunif(0,5)"), 
         detPriors=setPriors(intercept="dunif(-10, 10)", coefficient="dlogis(0, 1)", sd = "dunif(0,5)"),
         marginalized = FALSE, speciesID=quote(speciesID),
         modelInfo, .env){
 
    site_dim <- LHS[[3]] # sites, e.g. 1:N
    obs_dim <- LHS[[4]]  # obs, e.g. 1:J
    sp_dim <- LHS[[5]]   # species, e.g. 1:S

    y <- eval(LHS[[2]], envir=modelInfo$constants)
    S <- dim(y)[3]
    modelInfo$constants[[safeDeparse(speciesID)]] <- factor(1:S, levels=1:S)

    sf <- multispeciesFormula(stateformula, sp_dim, speciesID)
    df <- multispeciesFormula(detformula, sp_dim, speciesID)

    # Code for calculating occupancy and detection parameters
    param_code <- substitute({
      psi[SITEDIM, SPDIM] <- nimbleMacros::linPred(STATEFORM, link=logit, coefPrefix=STATEPREFIX, 
                                     sdPrefix=STATEPREFIX, priorSettings=STATEPRIOR, 
                                     center=SPID)  
      p[SITEDIM, OBSDIM, SPDIM] <- nimbleMacros::linPred(DETFORM, link=logit, coefPrefix=DETPREFIX, 
                                           sdPrefix=DETPREFIX, priorSettings=DETPRIOR, 
                                           center=SPID)
    }, 
    list(SITEDIM=site_dim, OBSDIM=obs_dim, SPDIM=sp_dim,
         STATEFORM=sf, DETFORM = df,
         STATEPREFIX=statePrefix, DETPREFIX = detPrefix,
         STATEPRIOR=statePriors, DETPRIOR=detPriors,
         SPID=speciesID)
    )
 
    if(!marginalized){
      # Latent state model    
      ymod_code <- substitute({
        z[SITEDIM, SPDIM] ~ nimbleMacros::forLoop(dbern(psi[SITEDIM, SPDIM]))

        RESP ~ nimbleMacros::forLoop(dbern(p[SITEDIM, OBSDIM, SPDIM]*z[SITEDIM, SPDIM]))
        }, 
        list(RESP=LHS, SITEDIM=site_dim, OBSDIM=obs_dim, SPDIM=sp_dim)
      )

      # Create initial values for z
      z <- apply(y, c(1,3), max, na.rm=TRUE)
      #z[z == 0] <- NA
      #modelInfo$constants$z <- z
      modelInfo$inits <- list(z = z)

    } else {
      resp <- LHS[[2]]
      site_idx <- as.name(modelInfo$indexCreator())
      sp_idx <- as.name(modelInfo$indexCreator())

      # Obs range for a specific site (e.g. 1:J[i])
      obs_dim2 <- replaceIndex(obs_dim, site_dim, site_idx) 
      nsamples <- obs_dim2[[3]]

      ymod_code <- substitute({
        for (SITEIDX in SITEDIM){
          for (SPIDX in SPDIM){
            RESP[SITEIDX, OBSDIM2, SPIDX] ~ dOcc_v(
              probOcc = psi[SITEIDX, SPIDX],
              probDetect = p[SITEIDX, OBSDIM2, SPIDX],
              len = NSAMP
              )
          }
        }

        }, 
        list(RESP=LHS[[2]], SITEIDX = site_idx, SPIDX = sp_idx, SITEDIM=site_dim, 
             OBSDIM=obs_dim, OBSDIM2 = obs_dim2, NSAMP=nsamples, SPDIM=sp_dim)
      )

    }

    code <- embedLinesInCurlyBrackets(list(param_code, ymod_code))

    # Return code and model info
    list(code=code, modelInfo=modelInfo)
  },
use3pieces=TRUE,
unpackArgs=TRUE
)
class(multispeciesOccupancy) <- "model_macro"

multispeciesFormula <- function(form, sp_idx, sp_vec){
  bars <- "||"
  if(form == ~1) bars <- "|"
  sp_vec <- safeDeparse(sp_vec)
  rand <- paste0("+ (", safeDeparse(form[[2]]), " ", bars, " ", 
                 sp_vec,"[",safeDeparse(sp_idx),"])")
  stats::as.formula(str2lang(paste(safeDeparse(form), rand)))
}

# not the same as nimble's version
safeDeparse <- function(inp) {
  out <- deparse(inp)
  paste(sapply(out, trimws), collapse=" ")
}
