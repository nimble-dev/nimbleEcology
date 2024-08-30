#' Internal helper nimbleFunctions for dNmixture distributions
#'
#' None of these functions should be called directly.
#'
#' @name dNmixture_steps
#'
#' @aliases dNmixture_BNB_steps dNmixture_BBP_steps dNmixture_BBNB_steps
#'   nimNmixPois_logFac
#'
#' @param x x from dNmixture distributions
#' @param lambda lambda from dNmixture distributions
#' @param theta theta from relevant dNmixture distributions
#' @param s s from relevant dNmixture distributions
#' @param Nmin start of summation over N
#' @param Nmax end of summation over N
#' @param sum_log_one_m_prob sum(log(1-prob)) from relevant dNmixture cases
#' @param sum_log_dbinom sum(log(dbinom(...))) from relevant dNmixture cases
#' @param sum_log_dbetabinom sum(log(dBetaBinom(...))) from relevant dNmixture
#'   cases
#' @param beta_m_x beta-x from relevant dNmixture cases
#' @param usingAD TRUE if called from one of the dNmixtureAD distributions
#' @param numN number of indices in the truncated sum for the N-mixture.
#' @param ff a derived vector of units calculated partway through the fast
#'   N-mixture algorithm.
#' @param max_index possibly the index of the max contribution to the summation.
#'   For AD cases this is set by heuristic. For non-AD cases it is -1 and will
#'   be determined automatically.
#' @details These are helper functions for the N-mixture calculations. They
#'   don't have an interpretation outside of that context and are not intended
#'   to be called directly.
#' @seealso \code{\link{dNmixture}}

#' @rdname dNmixture_steps
#' @export
nimNmixPois_logFac <- nimbleFunction(
  run = function(numN = integer(0),
                 ff = double(1),
                 max_index = integer(0, default=-1)) {
    fixed_max_index <- max_index != -1

    i <- 1L
    sum_ff_g1 <- 0
    if(!fixed_max_index) {
      if(numN == 1) {
        sum_ff_g1 <- ff[1]
        max_index <- 1
      } else {
        hit_pos <- FALSE
        while(i < numN & (ff[i] > 0 | !hit_pos)) {
          sum_ff_g1 <- sum_ff_g1 + ff[i]
          i <- i+1
          if (ff[i] > 0) {
            hit_pos <- TRUE
          }
        }
        max_index <- i-1
        if (ff[i] > 0 & numN != max_index + 1) {
          max_index <- i
          sum_ff_g1 <- sum_ff_g1 + ff[i]
        }
        if(max_index == 0 | !hit_pos) {
          max_index <- 1 # not sure this is relevant. it's defensive.
          sum_ff_g1 <- ff[1]
        }
      } #end numN != 1
    } else { # end !fixed_max_index
      # here we are in the case that max_index was provided, which we be in an AD case
      for(i in 1:max_index) {
        sum_ff_g1 <- sum_ff_g1 + ff[i]
      }
    }
    terms <- numeric(numN + 1)
    terms[max_index + 1] <- 1

    sumff <- sum_ff_g1 ## should be the same as sum(ff[1:max_index])

    for (i in 1:max_index) {
                                        # terms[i] <- 1 / exp(sum(ff[i:max_index]))
      terms[i] <- 1 / exp(sumff)
      sumff <- sumff - ff[i]
    }
    if(numN > 1) {
      sumff <- 0
      for (i in (max_index + 1):numN) {
                                        # terms[i + 1] <- exp(sum(ff[(max_index + 1):i]))
        sumff <- sumff + ff[i]
        terms[i + 1] <- exp(sumff)
      }
    }
    log_fac <- sum_ff_g1 + log(sum(terms)) # Final factor is the largest term * (all factors / largest term)    }
    return(log_fac)
    returnType(double())
  }, buildDerivs = list(run = list()))

#' @rdname dNmixture_steps
#' @importFrom stats dpois
#' @export
dNmixture_steps <- nimbleFunction(
  run = function(x = double(1),
                 lambda = double(),
                 Nmin = double(),
                 Nmax = double(),
                 sum_log_one_m_prob = double(),
                 sum_log_dbinom = double(),
                 usingAD = logical(0, default=FALSE)) {
    logProb <- -Inf
    max_index <- -1L
    NminC <- NmaxC <- 0L
    NminC <- ADbreak(Nmin)
    NmaxC <- ADbreak(Nmax)
    logProb <- dpois(Nmin, lambda, log = TRUE) + sum_log_dbinom
    if(Nmax > Nmin) {
      if(usingAD) {
        # We need a choice for the max_index in nimNmixPois_logFac that is not a function of parameters,
        # because it will be baked into the AD tape upon first call.
        # A heuristic guess is either 2*Nmin or 0.2 between Nmin and Nmax
        # whichever is smaller.
        # But actually summation starts at Nmin+1, so these are tweaked accordingly.
        # Note that mathematically the result should work for any choice
        # of max_index, and a good choice only provides
        # some stability against underflows.

        # for a couple of steps, max_index is on N scale, then re-used relative to Nmin
        max_index <- ADbreak(min(2*NminC,
                                 floor(NminC + 0.2*(NmaxC-NminC))))
        # Make completely sure we are at least one below the max, mean 2 below the max at this step
        max_index <- ADbreak(min(max_index,
                                 NmaxC - 2))
        # shift max_index relative to Nmin and 1-indexed
        max_index <- max_index - NminC + 1
        # And make completely sure we are not at 0. Not sure that could happen, but being defensive.
        if(max_index < 1) max_index <- 1
      }
      numN <- 0L
      numN <- NmaxC - NminC + 1 - 1  ## remember: +1 for the count, but -1 because the summation should run from N = maxN to N = minN + 1
      prods <- rep(0, numN)
      for (i in (NminC + 1):NmaxC) {
        prodi <- 1
        for (j in 1:length(x)){
          xj <- ADbreak(x[j])
          if(!is.na(xj)){
            prodi <- prodi * (i / (i - x[j]))
          }
        }
        prods[i - NminC] <- prodi / i
      }
      ff <- log(lambda) + sum_log_one_m_prob + log(prods)
      log_fac <- nimNmixPois_logFac(numN, ff, max_index)
      logProb <- logProb + log_fac
    }
    return(logProb)
    returnType(double())
  },
  buildDerivs = list(run = list(ignore = c("i","j","xj")))
)

##### N-mixture extensions #####
##### dNmixture_BNB_v #####
NULL
#' @rdname dNmixture_steps
#' @importFrom stats dnbinom
#' @export
dNmixture_BNB_steps <- nimbleFunction(
  run = function(x = double(1),
                 lambda = double(),
                 theta = double(),
                 Nmin = double(),
                 Nmax = double(),
                 sum_log_one_m_prob = double(),
                 sum_log_dbinom = double(),
                 usingAD = logical(0, default=FALSE)) {
    r <- 1 / theta
    pNB <- 1 / (1 + theta * lambda)
    logProb <- -Inf
    max_index <- -1L
    NminC <- NmaxC <- 0L
    NminC <- ADbreak(Nmin)
    NmaxC <- ADbreak(Nmax)
    logProb <- dnbinom(Nmin, size = r, prob = pNB, log = TRUE) + sum_log_dbinom
    if(Nmax > Nmin) {
      if(usingAD) {
        # see comments in basic case above
        max_index <- ADbreak(min(2*NminC,
                                 floor(NminC + 0.2*(NmaxC-NminC))))
        max_index <- ADbreak(min(max_index,
                                 NmaxC - 2))
        max_index <- max_index - NminC + 1
        if(max_index < 1) max_index <- 1
      }
      numN <- 0L
      numN <- NmaxC - NminC + 1 - 1 # remember...
      prods <- rep(0, numN)
      for (i in (NminC + 1):NmaxC){
        prods[i - NminC] <- (i + r - 1) * prod(i/(i - x)) / i

        prodi <- i + r - 1
        for (j in 1:length(x)){
          xj <- ADbreak(x[j])
          if(!is.na(xj)){
            prodi <- prodi * (i / (i - x[j]))
          }
        }
        prods[i - NminC] <- prodi / i
      }
      ff <- log(1 - pNB) + sum_log_one_m_prob + log(prods)
      log_fac <- nimNmixPois_logFac(numN, ff, max_index)
      logProb <- logProb + log_fac
    }
    return(logProb)
    returnType(double())
  },
  buildDerivs = list(run = list(ignore = c("i", "j", "xj")))
)


##### dNmixture_BBP_v #####
NULL
#' @rdname dNmixture_steps
#' @importFrom stats dpois
#' @export
dNmixture_BBP_steps <- nimbleFunction(
  run = function(x = double(1),
                 beta_m_x = double(1),
                 lambda = double(),
                 s = double(),
                 Nmin = double(),
                 Nmax = double(),
                 sum_log_dbetabinom = double(),
                 usingAD = logical(0, default=FALSE)) {
    logProb <- -Inf
    max_index <- -1L
    NminC <- NmaxC <- 0L
    NminC <- ADbreak(Nmin)
    NmaxC <- ADbreak(Nmax)
    logProb <- dpois(Nmin, lambda, log = TRUE) + sum_log_dbetabinom
    if(Nmax > Nmin) {
      if(usingAD) {
        # see comments in basic case above
        max_index <- ADbreak(min(2*NminC,
                                 floor(NminC + 0.2*(NmaxC-NminC))))
        max_index <- ADbreak(min(max_index,
                                 NmaxC - 2))
        max_index <- max_index - NminC + 1
        if(max_index < 1) max_index <- 1
      }
      numN <- 0L
      numN <- NmaxC - NminC + 1 - 1 # remember...
      prods <- rep(0, numN)
      # N.B. alpha+beta == s
      for (i in (NminC + 1):NmaxC){
        prodi <- 1
        for (j in 1:length(x)){
          xj <- ADbreak(x[j])
          if(!is.na(xj)){
            prodi <- prodi * (i * (i - 1 + beta_m_x[j]) / ((i - x[j]) * (s + i - 1)))
          }
          prods[i - NminC] <- prodi * (lambda / i)
        }
      }
      ff <- log(prods)
      log_fac <- nimNmixPois_logFac(numN, ff, max_index)
      logProb <- logProb + log_fac
    }
    return(logProb)
    returnType(double())
  },
  buildDerivs = list(run = list(ignore = c("i","j","xj")))
)


##### dNmixture_BBNB_v #####
NULL
#' @rdname dNmixture_steps
#' @importFrom stats dnbinom
#' @export
dNmixture_BBNB_steps <- nimbleFunction(
  run = function(x = double(1),
                 beta_m_x = double(1),
                 lambda = double(),
                 theta = double(),
                 s = double(),
                 Nmin = double(),
                 Nmax = double(),
                 sum_log_dbetabinom = double(),
                 usingAD = logical(0, default=FALSE)) {
    r <- 1 / theta
    pNB <- 1 / (1 + theta * lambda)
    logProb <- -Inf
    max_index <- -1L
    NminC <- NmaxC <- 0L
    NminC <- ADbreak(Nmin)
    NmaxC <- ADbreak(Nmax)
    logProb <- dnbinom(Nmin, size = r, prob = pNB, log = TRUE) + sum_log_dbetabinom
    if(Nmax > Nmin) {
      if(usingAD) {
        # see comments in basic case above
        max_index <- ADbreak(min(2*NminC,
                                 floor(NminC + 0.2*(NmaxC-NminC))))
        max_index <- ADbreak(min(max_index,
                                 NmaxC - 2))
        max_index <- max_index - NminC + 1
        if(max_index < 1) max_index <- 1
      }
      numN <- 0L
      numN <- NmaxC - NminC + 1 - 1 # remember...
      prods <- rep(0, numN)
        # N.B. alpha+beta == s
      for (i in (NminC + 1):NmaxC)
                prods[i - NminC] <- prod(i * (i - 1 + beta_m_x) / ((i - x) * (s + i - 1))) *
          ((1 - pNB) * (i + r - 1) / i)
      ff <- log(prods)
      log_fac <- nimNmixPois_logFac(numN, ff, max_index)
      logProb <- logProb + log_fac
    }
    return(logProb)
    returnType(double())
  },
  buildDerivs = list(run = list(ignore = c("i")))
)
