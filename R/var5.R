source("settings.R")
source("eqocean.R")
source("run.R")
source("plot.R")

q0 <- cos(k0 * x)
q2 <- cos(k2 * x)
q4 <- cos(k4 * x)
f <- matrix(0, imax, tmax)
f[, 1] <- 0.1 * cos(kf * x)

true_state <- run_forward(q0, q2, q4, f, tmax)

thobs <- 1:tmax
nthobs <- length(thobs)
sobs <- 0
htrue <- calc_h(true_state$q0, true_state$q2, true_state$q4)
hobs <- htrue[, thobs] + matrix(rnorm(imax * nthobs, 0, sobs), imax, nthobs)

tfobs <- c(1, 32) + 1
ntfobs <- length(tfobs)
sstr <- 0.01
fobs <- true_state$f[, tfobs] + matrix(rnorm(imax * ntfobs, 0, sstr), imax, ntfobs)

niter <- 10000
ftol <- 1e-5
ctol <- 1e-8
gtol <- 1e-8

smod <- 0.01
q0 <- sin(k0 * x)
q2 <- sin(k2 * x)
q4 <- sin(k4 * x)
sh <- 1
ss <- 1
lr <- 5e-3
ocost <- 1e8
ognorm <- 1
fcg = 1 # set to 1 to use CG
d <- numeric(3 * imax + imax * tmax)

dh <- matrix(0, imax, tmax)
ds <- matrix(0, imax, tmax)
f[, ] <- 0
for (ep in 1:niter) {
  state <- run_forward(q0, q2, q4, f, tmax, FALSE)
  h <- calc_h(state$q0, state$q2, state$q4)
  f <- state$f
  dh[, -thobs] <- 0
  dh[, thobs] <- hobs - h[, thobs]
  ds[, -tfobs] <- 0
  ds[, tfobs] <- fobs - f[, tfobs]
  adj <- run_adjoint(dh, ds)
  cost <- calc_cost(dh, ds, sh, ss)
  xf <- c(q0, q2, q4, f)
  g <- c(-adj$p0, -adj$p2, -adj$p4, adj$ps)
  gnorm <- sum(g^2)
  bta <- gnorm / ognorm
  d <- -g + fcg * bta * d
  xa <- xf + lr * d
  q0 <- xa[1:imax]
  q2 <- xa[1:imax + imax]
  q4 <- xa[1:imax + 2 * imax]
  f <- matrix(xa[(3 * imax + 1):length(xa)], imax, tmax)
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
state <- run_forward(q0, q2, q4, f, tmax, FALSE)
print_mse(state, true_state)
#saveRDS(state, "analysis_var5.rds")
plot_waves(x, state, true_state)
