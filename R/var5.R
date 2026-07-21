source("settings.R")
source("eqocean.R")
source("run.R")
source("plot.R")

q0 <- cos(k0 * x)
q2 <- cos(k2 * x)
q4 <- cos(k4 * x)
f <- 0.1 * cos(kf * x)

true_state <- run_forward(q0, q2, q4, f, tmax)

#tobs <- c(1, 16, 32) + 1
tobs <- 1:tmax
ntobs <- length(tobs)
sobs <- 0
htrue <- calc_h(true_state$q0, true_state$q2, true_state$q4)
hobs <- htrue[, tobs] + matrix(rnorm(imax * ntobs, 0, sobs), imax, ntobs)

tstr <- c(1, 32) + 1
tstr <- 2:tmax
#ntstr <- length(tstr)
sstr <- 0.01
fobs <- true_state$f[, tstr] + matrix(rnorm(imax * ntstr, 0, sstr), imax, ntstr)

niter <- 10000
ftol <- 1e-5
ctol <- 1e-8
gtol <- 1e-8

smod <- 0.01
q0 <- cos(k0 * x) + rnorm(imax, 0, smod)
q2 <- cos(k2 * x) + rnorm(imax, 0, smod)
q4 <- cos(k4 * x) + rnorm(imax, 0, smod)
#q0 <- sin(k0 * x)
#q2 <- sin(k2 * x)
#q4 <- sin(k4 * x)
sh <- 0.1
ss <- 1
lr <- 1e-2
decay_rate <- 1e-2
ocost <- 1e8
ognorm <- 1
fcg = 1 # set to 1 to use CG
d <- numeric(3 * imax + imax * tmax)

dh <- matrix(0, imax, tmax)
ds <- matrix(0, imax, tmax)
#fa <- matrix(0, imax, tmax)
fa <- matrix(rep(0.1 * cos(kf * x) + rnorm(imax, 0, smod)), imax, tmax)
#fa <- matrix(rep(0.1 * sin(kf * x), tmax), imax, tmax)
for (ep in 1:niter) {
  state <- run_forward(q0, q2, q4, fa, tmax)
  h <- calc_h(state$q0, state$q2, state$q4)
  fa <- state$f
  dh[, !tobs] <- 0
  dh[, tobs] <- hobs - h[, tobs]
  ds[, !tstr] <- 0
  ds[, tstr] <- fobs - state$f[, tstr]
  adj <- run_adjoint(dh, ds)
  cost <- calc_cost(dh, ds, sh, ss)
  xf <- c(q0, q2, q4, fa)
  g <- c(adj$p0, adj$p2, adj$p4, adj$ps)
  gnorm <- sum(g^2)
  bta <- gnorm / ognorm
  d <- -g + fcg * bta * d
#  if (ep > 1) lr <- lr * sum(dold * gold) / sum(g * d)
#  if (ep > 1) lr <- lr * sum((xf - xold)^2) / sum((xf - xold) * (g - gold))
#  if (ep > 1) lr <- lr * sum((xf - xold) * (g - gold)) / sum((g - gold)^2)
  lr <- lr / (1 + decay_rate * ep)
  xa <- xf + lr * d
#  xa <- si * xa + mu
  q0 <- xa[1:imax]
  q2 <- xa[1:imax + imax]
  q4 <- xa[1:imax + 2 * imax]
  fa <- matrix(xa[(3 * imax + 1):length(xa)], imax, tmax)
  f <- fa[, 1]
  dcost <- abs(cost - ocost)
  if ((ep - 1) %% 10 == 0) print_diag(ep, cost, gnorm)
  if (cost < ftol) {print_diag(ep, cost, gnorm, "stopping due to ftol"); break}
  if (dcost < ctol) {print_diag(ep, cost, gnorm, "stopping due to ctol"); break}
  if (gnorm < gtol) {print_diag(ep, cost, gnorm, "stopping due to gtol"); break}
  ocost <- cost
  ognorm <- gnorm
  gold <- g
  dold <- d
  xold <- xf
}
state <- run_forward(q0, q2, q4, f, tmax)
print_mse(state, true_state)
#saveRDS(state, "analysis_var5.rds")
plot_waves(x, state)
