
x <- runif(100, 0, 1)
c <- as.factor(rep(1:20, 5)); c_ran <- rnorm(20, 0, 0.5)
y <- rpois(100, exp(x * 3 + c_ran[c] + rnorm(100, 0, 0.1)))
practice_fit <- glmer(y ~ x + (1|c),  family = poisson())
class(practice_fit)
summary(practice_fit)

fit_Occ <- function(expr, data, method = c("MLE", "MCMC")) {




}
