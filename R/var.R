source("settings.R")
source("eqocean.R")
source("run.R")
source("plot.R")

q0 <- cos(k0 * x)
q2 <- cos(k2 * x)
q4 <- cos(k4 * x)
f <- 0.1 * cos(kf * x)

true_state <- run_forward(q0, q2, q4, f, tmax)

tobs <- c(0, 16, 32) + 1
ntobs <- length(tobs)
sobs <- 0.01
htrue <- calc_h(true_state$q0, true_state$q2, true_state$q4)
hobs <- htrue[, tobs] + matrix(rnorm(imax * ntobs, 0, sobs), imax, ntobs)

niter <- 1000
ftol <- 1e-5
ctol <- 1e-8
gtol <- 1e-8

smod <- 0.01
#q0 <- cos(k0 * x) + rnorm(imax, 0, smod)
#q2 <- cos(k2 * x) + rnorm(imax, 0, smod)
#q4 <- cos(k4 * x) + rnorm(imax, 0, smod)
q0 <- sin(k0 * x)
q2 <- sin(k2 * x)
q4 <- sin(k4 * x)
f <- 0.1 * cos(kf * x)
lr <- 0.1
ocost <- 1e8
ognorm <- 1
fcg = 1 # set to 1 to use CG
d <- numeric(3 * imax)

dh <- matrix(0, imax, tmax)
for (ep in 1:niter) {
  state <- run_forward(q0, q2, q4, f, tmax)
  h <- calc_h(state$q0, state$q2, state$q4)
  dh[, !tobs] <- 0
  dh[, tobs] <- h[, tobs] - hobs
  adj <- run_adjoint(dh)
  cost <- calc_cost(dh, 0)
  xf <- c(q0, q2, q4)
  g <- c(adj$p0, adj$p2, adj$p4)
  gnorm <- sum(g^2)
  bta <- gnorm / ognorm
  d <- -g + fcg * bta * d
  xa <- xf + lr * d
  q0 <- xa[1:imax]
  q2 <- xa[1:imax + imax]
  q4 <- xa[1:imax + 2 * imax]
  dcost <- abs(cost - ocost)
  if ((ep - 1) %% 10 == 0) print_diag(ep, cost, gnorm)
  if (cost < ftol) {print_diag(ep, cost, gnorm, "stopping due to ftol"); break}
  if (dcost < ctol) {print_diag(ep, cost, gnorm, "stopping due to ctol"); break}
  if (gnorm < gtol) {print_diag(ep, cost, gnorm, "stopping due to gtol"); break}
  ocost <- cost
  ognorm <- gnorm
}
state <- run_forward(q0, q2, q4, f, tmax)
#saveRDS(state, "analysis_var.rds")
plot_waves(x, state)
