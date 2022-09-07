
#' @rdname dNmixture_M
#'
#' @export
#'
dNmixture_MP_s <- nimbleFunction(
  run = function(x   = double(1),
                 mu  = double(),
                 p   = double(),
                 J   = double(),
                 log = integer(0, default = 0)) {

    x_tot <- sum(x)
    x_miss <- sum(x * seq(0, J - 1))


    # Mathematica version of MNB:
    ### loglik = Total[LogGamma[yrow+r]] - R*LogGamma[r] - ylogfact +
    ###   R*r*Log[r] + ytot*Log[mu] +
    ###   ytot*Log[p] + ysumj*Log[1-p] -
    ###   (ytot+R*r)*Log[r + mu*ptot];

    # Mathematica version of MP:
    ### loglik = ytot*Log[mu] +  ytot*Log[p]+ ysumj*Log[1-p] - R*mu*ptot-ylogfact;

    # term1   <- lgamma(r + x_tot) - lgamma(r) - sum(lfactorial(x))
    # term2   <- r * log(r) + x_tot * log(mu)
    # term3   <- x_tot * log(p) + x_miss * log(1 - p)
    # term4   <- -(x_tot + r) * log(r + mu * (1 - (1 - p) ^ J))
    logProb <- x_tot * log(mu) + x_tot * log(p) + x_miss * log(1 - p) -
      mu * (1 - (1 - p) ^ J) - sum(lfactorial(x))


    if (log) return(logProb)
    else return(exp(logProb))
    returnType(double())

  })

#' @rdname dNmixture_M
#'
#' @export
#'
rNmixture_MP_s <- nimbleFunction(
  run = function(n  = integer(),
                 mu = double(),
                 p  = double(),
                 J  = double()) {


    prob <- numeric(J + 1)
    for (i in 1:(J)) {
      prob[i] <- pow(1 - p, i - 1) * p
    }
    prob[J + 1] <- 1 - sum(prob[1:J])

    ans <- numeric(J + 1)
    n <- rpois(n = 1, lambda = mu)
    if (n > 0) {
      ans <- rmulti(n = 1, size = n, prob = prob)
    }

    return(ans[1:J])
    returnType(double(1))
  })


#' @rdname dNmixture_M
#'
#' @export
#'
dNmixture_MP_v <- nimbleFunction(
  run = function(x   = double(1),
                 mu  = double(),
                 p   = double(1),
                 J   = double(),
                 log = integer(0, default = 0)) {


    x_tot <- sum(x)
    x_miss <- sum(x * seq(0, J - 1))

    pp <- c(0, p)
    prob <- numeric(J)
    for (j in 1:J) {
      prob[j] <- prod(1 - pp[1:j]) * p[j]
    }
    ptot <- sum(prob)

    logProb <- x_tot * log(mu) + sum(x * log(prob)) -
      mu * ptot - sum(lfactorial(x))

    if (log) return(logProb)
    else return(exp(logProb))
    returnType(double())

  })


#' @rdname dNmixture_M
#'
#' @export
#'
rNmixture_MP_v <- nimbleFunction(
  run = function(n  = integer(),
                 mu = double(),
                 p  = double(1),
                 J  = double()) {


    pp <- c(0, p)
    prob <- numeric(J + 1)
    for (j in 1:J) {
      prob[j] <- prod(1 - pp[1:j]) * p[j]
    }
    prob[J + 1] <- 1 - sum(prob[1:J])

    ans <- numeric(J + 1)
    n <- rpois(n = 1, lambda = mu)
    if (n > 0) {
      ans <- rmulti(n = 1, size = n, prob = prob)
    }

    return(ans[1:J])
    returnType(double(1))
  })
