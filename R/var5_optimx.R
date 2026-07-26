library("optimx")
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
f[, ] <- 0
sh <- 1
ss <- 1

fn <- function(par, tmax, hobs, thobs, fobs, tfobs, sh, ss) {
  q0 <- par[1:imax]
  q2 <- par[1:imax + imax]
  q4 <- par[1:imax + 2 * imax]
  f <- matrix(par[(3 * imax + 1):length(par)], imax, tmax)
  state <- run_forward(q0, q2, q4, f, tmax, FALSE)
  h <- calc_h(state$q0, state$q2, state$q4)
  dh <- matrix(0, imax, tmax)
  dh[, thobs] <- hobs - h[, thobs]
  ds <- matrix(0, imax, tmax)
  ds[, tfobs] <- fobs - state$f[, tfobs]
  calc_cost(dh, ds, sh, ss)
}

gr <- function(par, tmax, hobs, thobs, fobs, tfobs, sh, ss) {
  q0 <- par[1:imax]
  q2 <- par[1:imax + imax]
  q4 <- par[1:imax + 2 * imax]
  f <- matrix(par[(3 * imax + 1):length(par)], imax, tmax)
  state <- run_forward(q0, q2, q4, f, tmax, FALSE)
  h <- calc_h(state$q0, state$q2, state$q4)
  dh <- matrix(0, imax, tmax)
  dh[, thobs] <- hobs - h[, thobs]
  ds <- matrix(0, imax, tmax)
  ds[, tfobs] <- fobs - state$f[, tfobs]
  adj <- run_adjoint(dh, ds)
  c(-adj$p0, -adj$p2, -adj$p4, adj$ps)
}

par <- c(q0, q2, q4, f)
alg <- "nvm"
cntl <- list(maxit = niter)
res <- optimr(par, method = alg, control = cntl,
              fn = function(par) fn(par, tmax, hobs, thobs, fobs, tfobs, sh, ss),
              gr = function(par) gr(par, tmax, hobs, thobs, fobs, tfobs, sh, ss))
q0 <- res$par[1:imax]
q2 <- res$par[1:imax + imax]
q4 <- res$par[1:imax + 2 * imax]
f <- matrix(res$par[(3 * imax + 1):length(res$par)], imax, tmax)
state <- run_forward(q0, q2, q4, f, tmax, FALSE)
print_mse(state, true_state)
#saveRDS(state, "analysis_var5.rds")
plot_waves(x, state, true_state)
