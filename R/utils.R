#' Helper function for fast N-mixture calculation
#' @name nimNmixPois_logFac
#' @param numN first argument for helper function nimNmixPois_logFac,
#'   representing the number of indices in the truncated sum for the N-mixture.
#' @param ff second argument for helper function nimNmixPois_logFac, a derived
#'   vector of units calculated partway through the fast N-mixture algorithm.
#' @rdname nimNmixPois_logFac
#' @details This is a helper function for the fast N-mixture calculation. It
#'   runs an iterative calculation present in all N-mixture varieties. It
#'   doesn't have an interpretation outside of that context.
#' @export
#' @seealso \code{\link{dNmixture}}
##### nimNmixPois_logFac #####
nimNmixPois_logFac <- nimbleFunction(
  run = function(numN = double(0),
                 ff = double(1)) {
    i <- 1
    sum_ff_g1 <- 0
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
    
    ## terms <- numeric(numN + 1) ## This cannot compile with AD
    terms_len <- ADbreak(numN) + 1
    terms <- numeric(length = ADbreak(terms_len), value = 0)
    terms[max_index + 1] <- 1
    
    sumff <- sum_ff_g1 ## should be the same as sum(ff[1:max_index])
    
    for (i in 1:max_index) {
      # terms[i] <- 1 / exp(sum(ff[i:max_index]))
      terms[i] <- 1 / exp(sumff)
      sumff <- sumff - ff[i]
    }
    
    sumff <- 0
    for (i in (max_index + 1):ADbreak(numN)) {
      # terms[i + 1] <- exp(sum(ff[(max_index + 1):i]))
      sumff <- sumff + ff[i]
      terms[i + 1] <- exp(sumff)
    }
    log_fac <- sum_ff_g1 + log(sum(terms)) # Final factor is the largest term * (all factors / largest term)    }
    return(log_fac)
    returnType(double())
  },
  buildDerivs = list(run = list(ignore = c("i", "max_index")))
)
