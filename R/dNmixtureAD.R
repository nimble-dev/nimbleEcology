# dNmixtureAD
#' N-mixture distributions with AD support for use in \code{nimble} models
#'
#' \code{dNmixtureAD_s} and \code{dNmixtureAD_v} provide Poisson-Binomial
#' mixture distributions of abundance ("N-mixture") for use in \code{nimble}
#' models when automatic differentiation may be needed by an algorithm.
#' Overdispersion alternatives are also provided.
#'
#' @name dNmixtureAD
#' @aliases dNmixtureAD_s dNmixtureAD_v rNmixtureAD_s rNmixtureAD_v dNmixtureAD_BNB_v
#'   dNmixtureAD_BNB_s dNmixtureAD_BNB_oneObs dNmixtureAD_BBP_v dNmixtureAD_BBP_s
#'   dNmixtureAD_BBP_oneObs dNmixtureAD_BBNB_v dNmixtureAD_BBNB_s
#'   rNmixtureAD_BBNB_oneObs rNmixtureAD_BNB_v rNmixtureAD_BNB_s rNmixtureAD_BNB_oneObs
#'   rNmixtureAD_BBP_v rNmixtureAD_BBP_s rNmixtureAD_BBP_oneObs rNmixtureAD_BBNB_v
#'   rNmixtureAD_BBNB_s rNmixtureAD_BBNB_oneObs
#'
#' @author Ben Goldstein, Lauren Ponisio, and Perry de Valpine
#'
#' @param x vector of integer counts from a series of sampling occasions.
#' @param lambda expected value of the Poisson distribution of true abundance
#' @param theta abundance overdispersion parameter required for negative binomial
#'        (*NB) N-mixture models. theta is parameterized such that variance of
#'        the negative binomial variable x is \code{lambda^2 * theta + lambda}
#' @param prob detection probability (scalar for \code{dNmixture_s}, vector for \code{dNmixture_v}).
#' @param s detection overdispersion parameter required for beta binomial (BB*)
#'        N-mixture models. s is parameterized such that variance of the beta
#'        binomial variable x is \code{V(x) = N \* prob \* (1-prob) \* (N +
#'        s) / (s + 1)}
#' @param Nmin minimum abundance to sum over for the mixture probability. Must be provided.
#' @param Nmax maximum abundance to sum over for the mixture probability. Must be provided.
#' @param len The length of the x vector
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return probability.
#' @param n number of random draws, each returning a vector of length
#'        \code{len}. Currently only \code{n = 1} is supported, but the
#'        argument exists for standardization of "\code{r}" functions.
#'
#' @details These nimbleFunctions provide distributions that can be
#'     used directly in R or in \code{nimble} hierarchical models (via
#'     \code{\link[nimble]{nimbleCode}} and
#'     \code{\link[nimble]{nimbleModel}}).
#'
#' See \code{\link{dNmixture}} for more information about the N-mixture
#' distributions.
#'
#' The versions here can be used in models that will be used by algorithms that
#' use nimble's system for automatic differentiation (AD). The primary
#' difference is that \code{Nmin} and \code{Nmax} must be provided. There are no
#' automatic defaults for these.
#'
#' In the AD system some kinds of values are "baked in" (cannot be changed) to
#' the AD calculations from the first call, unless and until the AD calculations
#' are reset. For all variants of the \code{dNmixtureAD} distributions, the
#' sizes of the inputs as well as \code{Nmin} and \code{Nmax} are baked in.
#' These can be different for different iterations through a for loop (or nimble
#' model declarations with different indices, for example), but the sizes and
#' \code{Nmin} and \code{Nmax} values for each specific iteration will be
#' "baked in" after the first call.
#'
#' @return The probability (or likelihood) or log probability of an observation
#'   vector \code{x}.

##### Regular N-mixture #####
NULL
#' @rdname dNmixtureAD
#' @export
dNmixtureAD_v <- nimbleFunction(
  run = function(x = double(1),
                 lambda = double(),
                 prob = double(1),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double(),
                 log = integer(0, default = 0)) {
    if (length(x) != len) stop("in dNmixtureAD_v, len must equal length(x).")
    if (len != length(prob)) stop("in dNmixtureAD_v, len must equal length(prob).")
    if (lambda < 0)
      if (log) return(-Inf) else return(0)
    if ((Nmin == -1) | (Nmax == -1))
      stop("Must provide Nmin and Nmax in AD version of dNmixture distributions")

    max_x <- 0
    for (i in 1:length(x)){
      xi <- ADbreak(x[i])
      if(!is.na(xi)){
        if(x[i] > max_x) max_x <- x[i]
      }
    }
    Nmin <- ADbreak(max( max_x, Nmin ))

    sum_log_dbinom <- 0
    sum_log_one_m_prob <- 0
    any_not_na <- FALSE
    # You can't combine this with the loop above because we need to know the final Nmin
    for (i in 1:length(x)){
      xi <- ADbreak(x[i])
      if(!is.na(xi)){
        sum_log_one_m_prob <- sum_log_one_m_prob + log(1 - prob[i])
        sum_log_dbinom <- sum_log_dbinom + dbinom(x[i], size = Nmin, prob = prob[i], log=TRUE)
        any_not_na <- TRUE
      }
    }

    logProb <- 0
    if(any_not_na){
      logProb <- dNmixture_steps(x, lambda, Nmin, Nmax, sum_log_one_m_prob,
                                 sum_log_dbinom, usingAD=TRUE)
    }
    if (log) return(logProb)
    else return(exp(logProb))
    returnType(double())
  }, buildDerivs = list(run = list(ignore=c("i","xi")))
)

#' @rdname dNmixtureAD
#' @export
dNmixtureAD_s <- nimbleFunction(
  run = function(x = double(1),
                 lambda = double(),
                 prob = double(0),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double(),
                 log = integer(0, default = 0)) {
    if (length(x) != len) stop("in dNmixtureAD_v, len must equal length(x).")
    if (lambda < 0)
      if (log) return(-Inf) else return(0)
    if ((Nmin == -1) | (Nmax == -1))
      stop("Must provide Nmin and Nmax in AD version of dNmixture distributions")
    
    max_x <- 0
    for (i in 1:length(x)){
      xi <- ADbreak(x[i])
      if(!is.na(xi)){
        if(x[i] > max_x) max_x <- x[i]
      }
    }
    Nmin <- ADbreak(max( max_x, Nmin ))

    sum_log_dbinom <- 0
    sum_log_one_m_prob <- 0
    any_not_na <- FALSE
    for (i in 1:length(x)){
      xi <- ADbreak(x[i])
      if(!is.na(xi)){
        sum_log_one_m_prob <- sum_log_one_m_prob + log(1 - prob)
        sum_log_dbinom <- sum_log_dbinom + dbinom(x[i], size = Nmin, prob = prob, log=TRUE)
        any_not_na <- TRUE
      }
    }

    logProb <- 0
    if(any_not_na){
      logProb <- dNmixture_steps(x, lambda, Nmin, Nmax, sum_log_one_m_prob,
                                 sum_log_dbinom, usingAD=TRUE)
    }
    if (log) return(logProb)
    else return(exp(logProb))
    returnType(double())
  }, buildDerivs = list(run = list(ignore=c("i","xi")))
)

NULL
#' @rdname dNmixtureAD
#' @export
rNmixtureAD_v <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 prob = double(1),
                 Nmin = double(0),
                 Nmax = double(0),
                 len = double()) {
    if ((Nmin == -1) | (Nmax == -1))
      stop("Must provide Nmin and Nmax in AD version of dNmixture distributions")
    return(rNmixture_v(n,lambda,prob,Nmin,Nmax,len))
    returnType(double(1))
  })

#' @rdname dNmixtureAD
#' @export
rNmixtureAD_s <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 prob = double(0),
                 Nmin = double(0),
                 Nmax = double(0),
                 len = double()) {
    if ((Nmin == -1) | (Nmax == -1))
      stop("Must provide Nmin and Nmax in AD version of dNmixture distributions")
    return(rNmixture_s(n,lambda,prob,Nmin,Nmax,len))
    returnType(double(1))
  })

##### BNB cases #####
#' @rdname dNmixtureAD
#' @export
dNmixtureAD_BNB_v <- nimbleFunction(
  run = function(x = double(1),
                 lambda = double(),
                 theta = double(),
                 prob = double(1),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double(),
                 log = integer(0, default = 0)) {
    if (length(x) != len) stop("in dNmixtureAD_BNB_v, len must equal length(x).")
    if (len != length(prob)) stop("in dNmixtureAD_BNB_v, len must equal length(prob).")
    if (theta <= 0)
      if (log) return(-Inf) else return(0)
    if (lambda < 0)
      if (log) return(-Inf) else return(0)
    if ((Nmin == -1) | (Nmax == -1))
      stop("Must provide Nmin and Nmax in AD version of dNmixture distributions")
    
    max_x <- 0
    for (i in 1:length(x)){
      xi <- ADbreak(x[i])
      if(!is.na(xi)){
        if(x[i] > max_x) max_x <- x[i]
      }
    }

    Nmin <- ADbreak(max( max_x, Nmin ))

    sum_log_dbinom <- 0
    sum_log_one_m_prob <- 0
    any_not_na <- FALSE
    for (i in 1:length(x)){
      xi <- ADbreak(x[i])
      if(!is.na(xi)){
        sum_log_one_m_prob <- sum_log_one_m_prob + log(1 - prob[i])
        sum_log_dbinom <- sum_log_dbinom + dbinom(x[i], size = Nmin, prob = prob[i], log=TRUE)
        any_not_na <- TRUE
      }
    }

    logProb <- 0
    if(any_not_na){
      logProb <- dNmixture_BNB_steps(x,lambda,theta,Nmin,Nmax, sum_log_one_m_prob,
                                     sum_log_dbinom, usingAD = TRUE)
    }
    if (log) return(logProb)
    else return(exp(logProb))
    returnType(double())
  }, buildDerivs = list(run = list(ignore=c("i", "xi")))
)

#' @rdname dNmixtureAD
#' @export
dNmixtureAD_BNB_s <- nimbleFunction(
  run = function(x = double(1),
                 lambda = double(),
                 theta = double(),
                 prob = double(0),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double(),
                 log = integer(0, default = 0)) {
    if (length(x) != len) stop("in dNmixtureAD_BNB_v, len must equal length(x).")
    if (theta <= 0)
      if (log) return(-Inf) else return(0)
    if (lambda < 0)
      if (log) return(-Inf) else return(0)
    if ((Nmin == -1) | (Nmax == -1))
      stop("Must provide Nmin and Nmax in AD version of dNmixture distributions")
    
    max_x <- 0
    for (i in 1:length(x)){
      xi <- ADbreak(x[i])
      if(!is.na(xi)){
        if(x[i] > max_x) max_x <- x[i]
      }
    }

    Nmin <- ADbreak(max( max_x, Nmin ))

    sum_log_dbinom <- 0
    sum_log_one_m_prob <- 0
    any_not_na <- FALSE
    for (i in 1:length(x)){
      xi <- ADbreak(x[i])
      if(!is.na(xi)){
        sum_log_one_m_prob <- sum_log_one_m_prob + log(1 - prob)
        sum_log_dbinom <- sum_log_dbinom + dbinom(x[i], size = Nmin, prob = prob, log=TRUE)
        any_not_na <- TRUE
      }
    }

    logProb <- 0
    if(any_not_na){
      logProb <- dNmixture_BNB_steps(x,lambda,theta,Nmin,Nmax, sum_log_one_m_prob,
                                     sum_log_dbinom, usingAD = TRUE)
    }
    if (log) return(logProb)
    else return(exp(logProb))
    returnType(double())
  }, buildDerivs = list(run = list(ignore=c("i", "xi")))
)

#' @rdname dNmixtureAD
#' @export
dNmixtureAD_BNB_oneObs <- nimbleFunction(
  run = function(x = double(),
                 lambda = double(),
                 theta = double(),
                 prob = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 log = integer(0, default = 0)) {
    xvec <- numeric(value = x, length=1)
    return(dNmixtureAD_BNB_s(xvec,lambda,theta,prob,Nmin,Nmax,1,log))
    returnType(double())
  }, buildDerivs = list(run = list())
)

#' @rdname dNmixtureAD
#' @export
rNmixtureAD_BNB_oneObs <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 theta = double(),
                 prob = double(0),
                 Nmin = double(0),
                 Nmax = double(0)) {
    if ((Nmin == -1) | (Nmax == -1))
      stop("Must provide Nmin and Nmax in AD version of dNmixture distributions")
    return(rNmixture_BNB_oneObs(n,lambda,theta,prob,Nmin,Nmax))
    returnType(double(1))
  })

### BBP cases ###
#' @rdname dNmixtureAD
#' @export
dNmixtureAD_BBP_v <- nimbleFunction(
  run = function(x = double(1),
                 lambda = double(),
                 prob = double(1),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double(),
                 log = integer(0, default = 0)) {
    if (length(x) != len) stop("in dNmixtureAD_BBP_v, len must equal length(x).")
    if (len != length(prob)) stop("in dNmixtureAD_BBP_v, len must equal length(prob).")
    if (s <= 0)
      if (log) return(-Inf) else return(0)
    alpha <- prob * s
    beta <- s - prob * s
    if (lambda < 0)
      if (log) return(-Inf) else return(0)
    ## For beta-binomial N-mixtures , the conditional distribution of (N - x |
    ## x) doesn't have a nice closed-form expression.
    if (Nmin == -1 | Nmax == -1) {
      stop("Dynamic choice of Nmin/Nmax is not supported for beta binomial N-mixtures.")
    }
    max_x <- 0
    any_not_na <- FALSE
    for (i in 1:length(x)){
      xi <- ADbreak(x[i])
      if(!is.na(xi)){
        if(x[i] > max_x) max_x <- x[i]
        any_not_na <- TRUE
      }
    }
    Nmin <- ADbreak(max( max_x, Nmin )) ## set Nmin to at least the largest x

    logProb <- 0
    if(any_not_na){
      logProb <- dNmixture_BBP_steps(x, beta-x, lambda, s, Nmin, Nmax,
                                     dBetaBinom_v(x, Nmin, alpha, beta, len, log = TRUE), usingAD=TRUE)
    }
    if (log) return(logProb)
    else return(exp(logProb))
    returnType(double())
  }, buildDerivs = list(run = list(ignore=c("i", "xi")))
)

#' @rdname dNmixtureAD
#' @export
dNmixtureAD_BBP_s <- nimbleFunction(
  run = function(x = double(1),
                 lambda = double(),
                 prob = double(),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double(),
                 log = integer(0, default = 0)) {
    if (length(x) != len) stop("in dNmixtureAD_BBP_s, len must equal length(x).")
    if (s <= 0)
      if (log) return(-Inf) else return(0)
    alpha <- prob * s
    beta <- s - prob * s
    if (lambda < 0)
      if (log) return(-Inf) else return(0)
    if (Nmin == -1 | Nmax == -1) {
      stop("Dynamic choice of Nmin/Nmax is not supported for beta binomial N-mixtures.")
    }
    #Clen <- 0L
    #Clen <- ADbreak(len)
    max_x <- 0
    any_not_na <- FALSE
    for (i in 1:length(x)){
      xi <- ADbreak(x[i])
      if(!is.na(xi)){
        if(x[i] > max_x) max_x <- x[i]
        any_not_na <- TRUE
      }
    }
    Nmin <- ADbreak(max( max_x, Nmin )) ## set Nmin to at least the largest x
    
    logProb <- 0
    if(any_not_na){
      logProb <- dNmixture_BBP_steps(x, beta-x, lambda, s, Nmin, Nmax,
                                     dBetaBinom_s(x, Nmin, alpha, beta, len, log = TRUE), usingAD=TRUE)
    }
    if (log) return(logProb)
    else return(exp(logProb))
    returnType(double())
  }, buildDerivs = list(run=list(ignore=c("i", "xi")))
)

#' @rdname dNmixtureAD
#' @export
dNmixtureAD_BBP_oneObs <- nimbleFunction(
  run = function(x = double(),
                 lambda = double(),
                 prob = double(),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 log = integer(0, default = 0)) {
    xvec <- numeric(value = x, length = 1)
    return(dNmixtureAD_BBP_s(xvec,lambda,prob,s,Nmin,Nmax,1,log))
    returnType(double())
  }, buildDerivs=list(run=list())
)

## BBNB cases ##
#' @rdname dNmixtureAD
#' @export
dNmixtureAD_BBNB_v <- nimbleFunction(
  run = function(x = double(1),
                 lambda = double(),
                 theta = double(),
                 prob = double(1),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double(),
                 log = integer(0, default = 0)) {
    if (length(x) != len) stop("in dNmixtureAD_BBNB_v, len must equal length(x).")
    if (len != length(prob)) stop("in dNmixtureAD_BBNB_v, len must equal length(prob).")
    if (s <= 0)
      if (log) return(-Inf) else return(0)
    if (theta <= 0)
      if (log) return(-Inf) else return(0)
    alpha <- prob * s
    beta <- s - prob * s
    if (lambda < 0)
      if (log) return(-Inf) else return(0)
    ## see comments above
    if (Nmin == -1 | Nmax == -1) {
      stop("Dynamic choice of Nmin/Nmax is not supported for beta binomial N-mixtures.")
    }

    max_x <- 0
    any_not_na <- FALSE
    for (i in 1:length(x)){
      xi <- ADbreak(x[i])
      if(!is.na(xi)){
        if(x[i] > max_x) max_x <- x[i]
        any_not_na <- TRUE
      }
    }
    Nmin <- ADbreak(max( max_x, Nmin )) ## set Nmin to at least the largest x

    logProb <- 0
    if(any_not_na){
      logProb <- dNmixture_BBNB_steps(x, beta-x, lambda, theta, s, Nmin, Nmax,
                                     dBetaBinom_v(x, Nmin, alpha, beta, len, log = TRUE), usingAD=TRUE)
    }
    if (log) return(logProb)
    else return(exp(logProb))
    returnType(double())
  }, buildDerivs=list(run=list(ignore=c("i", "xi")))
)

#' @rdname dNmixtureAD
#' @export
dNmixtureAD_BBNB_s <- nimbleFunction(
  run = function(x = double(1),
                 lambda = double(),
                 theta = double(),
                 prob = double(),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double(),
                 log = integer(0, default = 0)) {
    if (length(x) != len) stop("in dNmixtureAD_BBNB_s, len must equal length(x).")
    if (s <= 0)
      if (log) return(-Inf) else return(0)
    if (theta <= 0)
      if (log) return(-Inf) else return(0)
    r <- 1 / theta
    pNB <- 1 / (1 + theta * lambda)
    alpha <- prob * s
    beta <- s - prob * s
    if (lambda < 0)
      if (log) return(-Inf) else return(0)
    ## See comments above
    if (Nmin == -1 | Nmax == -1) {
      stop("Dynamic choice of Nmin/Nmax is not supported for beta binomial N-mixtures.")
    }
#    Clen <- 0L
#    Clen <- ADbreak(len)
    max_x <- 0
    any_not_na <- FALSE
    for (i in 1:length(x)){
      xi <- ADbreak(x[i])
      if(!is.na(xi)){
        if(x[i] > max_x) max_x <- x[i]
        any_not_na <- TRUE
      }
    }
    Nmin <- ADbreak(max( max_x, Nmin )) ## set Nmin to at least the largest x
    
    logProb <- 0
    if(any_not_na){
      logProb <- dNmixture_BBNB_steps(x, beta-x, lambda, theta, s, Nmin, Nmax,
                                     dBetaBinom_s(x, Nmin, alpha, beta, len, log = TRUE), usingAD=TRUE)
    }
    if (log) return(logProb)
    else return(exp(logProb))
    returnType(double())
  }, buildDerivs=list(run=list(ignore=c("i", "xi")))
)

#' @rdname dNmixtureAD
#' @export
dNmixtureAD_BBNB_oneObs <- nimbleFunction(
  run = function(x = double(),
                 lambda = double(),
                 theta = double(),
                 prob = double(),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 log = integer(0, default = 0)) {
    xvec <- numeric(x, length = 1)
    return(dNmixtureAD_BBNB_s(xvec, lambda, theta, prob, s, Nmin, Nmax, 1, log))
    returnType(double())
  }, buildDerivs=list(run=list())
)

##### rNmixtureAD extensions #####
NULL
#' @rdname dNmixtureAD
#' @export
rNmixtureAD_BNB_v <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 theta = double(),
                 prob = double(1),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double()) {
    return(rNmixture_BNB_v(n,lambda,theta,prob,Nmin,Nmax,len))
    returnType(double(1))
  })

NULL
#' @rdname dNmixtureAD
#' @export
rNmixtureAD_BNB_s <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 theta = double(),
                 prob = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double()) {
    return(rNmixture_BNB_s(n,lambda,theta,prob,Nmin,Nmax,len))
    returnType(double(1))
  })


NULL
#' @rdname dNmixtureAD
#' @export
rNmixtureAD_BNB_oneObs <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 theta = double(),
                 prob = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1)) {
    return(rNmixture_BNB_oneObs(n,lambda,theta,prob,Nmin,Nmax))
    returnType(double())
  })

NULL
#' @rdname dNmixtureAD
#' @export
rNmixtureAD_BBP_v <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 prob = double(1),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double()) {
    return(rNmixture_BBP_v(n,lambda,prob,s,Nmin,Nmax,len))
    returnType(double(1))
  })

NULL
#' @rdname dNmixtureAD
#' @export
rNmixtureAD_BBP_s <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 prob = double(),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double()) {
    return(rNmixture_BBP_s(n,lambda,prob,s,Nmin,Nmax,len))
    returnType(double(1))
  })

NULL
#' @rdname dNmixtureAD
#' @export
rNmixtureAD_BBP_oneObs <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 prob = double(),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1)) {
    return(rNmixture_BBP_oneObs(n,lambda,prob,s,Nmin,Nmax))
    returnType(double())
  })

NULL
#' @rdname dNmixtureAD
#' @export
rNmixtureAD_BBNB_v <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 theta = double(),
                 prob = double(1),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double()) {
    return(rNmixture_BBNB_v(n,lambda,theta,prob,s,Nmin,Nmax,len))
    returnType(double(1))
  })

NULL
#' @rdname dNmixtureAD
#' @export
rNmixtureAD_BBNB_s <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 theta = double(),
                 prob = double(),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1),
                 len = double()) {
    return(rNmixture_BBNB_s(n,lambda,theta,prob,s,Nmin,Nmax,len))
    returnType(double(1))
  })

NULL
#' @rdname dNmixtureAD
#' @export
rNmixtureAD_BBNB_oneObs <- nimbleFunction(
  run = function(n = double(),
                 lambda = double(),
                 theta = double(),
                 prob = double(),
                 s = double(),
                 Nmin = double(0, default = -1),
                 Nmax = double(0, default = -1)) {
    return(rNmixture_BBNB_oneObs(n,lambda,theta,prob,s,Nmin,Nmax))
    returnType(double())
  })
