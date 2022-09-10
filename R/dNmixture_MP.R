
#' @rdname dNmixture_M
#'
#' @export
#'
dNmixture_MP_s <- nimbleFunction(
  run = function(x       = double(1),
                 lambda  = double(),
                 p       = double(),
                 J       = double(),
                 log     = integer(0, default = 0)) {

    x_tot <- sum(x)
    x_miss <- sum(x * seq(0, J - 1))

    p_vec <- rep(p, J)

    logProb <- x_tot * log(lambda) + sum(x * log(p_vec)) -
      lambda * sum(p_vec) - sum(lfactorial(x))


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
                 lambda = double(),
                 p  = double(),
                 J  = double()) {

    prob <- c(rep(p, J), 1 - (p * J))

    ans <- numeric(J + 1)
    n <- rpois(n = 1, lambda = lambda)
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
                 lambda  = double(),
                 p   = double(1),
                 J   = double(),
                 log = integer(0, default = 0)) {


    x_tot <- sum(x)
    x_miss <- sum(x * seq(0, J - 1))

    logProb <- x_tot * log(lambda) + sum(x * log(p)) -
      lambda * sum(p) - sum(lfactorial(x))

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
                 lambda = double(),
                 p  = double(1),
                 J  = double()) {

    prob <- c(p, 1 - sum(p))

    ans <- numeric(J + 1)
    n <- rpois(n = 1, lambda = lambda)
    if (n > 0) {
      ans <- rmulti(n = 1, size = n, prob = prob)
    }

    return(ans[1:J])
    returnType(double(1))
  })
