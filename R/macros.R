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
#' @author Ken Kellner
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
    
    # Code to calculate parameters psi and p
    param_code <- substitute({
      psi[SITEDIM] <- linPred(STATEFORM, link=logit, coefPrefix=STATEPREFIX, 
                                sdPrefix=STATEPREFIX, priorSettings=STATEPRIOR, 
                                centerVar=CENTVAR)  
      p[SITEDIM, OBSDIM] <- linPred(DETFORM, link=logit, coefPrefix=DETPREFIX, 
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
        z[SITEDIM] ~ forLoop(dbern(psi[SITEDIM]))
        Y ~ forLoop(dbern(p[SITEDIM, OBSDIM]*z[SITEDIM]))
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

      site_idx <- as.name(modelInfo$indexCreator())

      # Also check if integer
      obs_dim2 <- replaceIndex(obs_dim, site_dim, site_dix) 
      nsamples <- obs_dim2[[3]]

      ymod_code <- substitute({  
        for (SITEIDX in SITEDIM){
          Y[SITEIDX, OBSDIM2] ~ dOcc_v(
            probOcc = psi[SITEIDX],
            probDetect = p[SITEIDX, OBSDIM2],
            len = NSAMP
            )
        }
        }, 
        list(RESP=LHS[[2]], SITEDIM=site_dim, OBSDIM=obs_dim,
             SITEIDX=site_idx, OBSDIM2=obs_dim2, NSAMP=nsamples)
      )
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
  assign("forLoop", nimbleMacros::forLoop, envir = globalenv())
  assign("linPred", nimbleMacros::linPred, envir = globalenv())
  assign("priors", nimbleMacros::priors, envir = globalenv())
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
