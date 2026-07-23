library("optimx")

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
sobs <- 0
htrue <- calc_h(true_state$q0, true_state$q2, true_state$q4)
hobs <- htrue[, tobs] + matrix(rnorm(imax * ntobs, 0, sobs), imax, ntobs)

smod <- 0.01
#q0 <- cos(k0 * x) + rnorm(imax, 0, smod)
#q2 <- cos(k2 * x) + rnorm(imax, 0, smod)
#q4 <- cos(k4 * x) + rnorm(imax, 0, smod)
q0 <- sin(k0 * x)
q2 <- sin(k2 * x)
q4 <- sin(k4 * x)
f <- 0.1 * cos(kf * x)
lr <- 0.1

fn <- function(par, f, tmax, hobs, tobs) {
  q0 <- par[1:imax]
  q2 <- par[1:imax + imax]
  q4 <- par[1:imax + 2 * imax]
  state <- run_forward(q0, q2, q4, f, tmax)
  h <- calc_h(state$q0, state$q2, state$q4)
  dh <- matrix(0, imax, tmax)
  dh[, tobs] <- hobs - h[, tobs]
  calc_cost(dh, 0)
}

gr <- function(par, f, tmax, hobs, tobs) {
  q0 <- par[1:imax]
  q2 <- par[1:imax + imax]
  q4 <- par[1:imax + 2 * imax]
  state <- run_forward(q0, q2, q4, f, tmax)
  h <- calc_h(state$q0, state$q2, state$q4)
  dh <- matrix(0, imax, tmax)
  dh[, tobs] <- hobs - h[, tobs]
  adj <- run_adjoint(dh)
  c(adj$p0, adj$p2, adj$p4)
}

par <- c(q0, q2, q4)
alg <- "nvm"
cntl <- list(maxit = 1000, tol = 1e-5)
res <- optimr(par, method = alg, control = cntl,
              fn = function(par) fn(par, f, tmax, hobs, tobs),
              gr = function(par) gr(par, f, tmax, hobs, tobs))
q0 <- res$par[1:imax]
q2 <- res$par[1:imax + imax]
q4 <- res$par[1:imax + 2 * imax]
state <- run_forward(q0, q2, q4, f, tmax)
#saveRDS(state, "analysis_var.rds")
plot_waves(x, state, true_state)
