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

s_obs <- 0.01
rmat <- diag(imax) / sobs^2
hmat <- matrix(c(rep(1, imax), numeric(imax * (imax - 1)),
                 rep(1, imax), numeric(imax * (imax - 1)),
                 rep(-1, imax), numeric(imax * (imax - 1))), imax, 3 * imax)

nmem <- 10

s_mod <- 0.01
s_ens <- 0.1
#q0 <- cos(k0 * x) + rnorm(imax, 0, s_mod)
#q2 <- cos(k2 * x) + rnorm(imax, 0, s_mod)
#q4 <- cos(k4 * x) + rnorm(imax, 0, s_mod)
q0 <- sin(k0 * x)
q2 <- sin(k2 * x)
q4 <- sin(k4 * x)
xf <- c(q0, q2, q4)
xe <- matrix(rnorm(3 * imax * nmem, 0, s_ens), 3 * imax, nmem)
pmat <- xe %*% t(xe) / (nmem - 1)
f <- 0.1 * cos(kf * x)

ensemble_update <- function(xf, xe, yo, hmat, rmat, infl = 0.1, std = 0.01) {
  nmem <- ncol(xe)
  omax <- length(yo)
  xe <- sqrt(1 + infl) * xe
  pmat <- xe %*% t(xe) / (nmem - 1)
  kmat <- pmat %*% t(hmat) %*% solve(hmat %*% pmat %*% t(hmat) + rmat)
  xf <- xf + kmat %*% (yo - hmat %*% xf)
  ye <- matrix(rnorm(omax * nmem, 0, std), omax, nmem)
  xe <- xe + kmat %*% (ye - hmat %*% xe)
  print(sum(diag(pmat)))
  list(xf = xf, xe = xe)
}

dtobs <- diff(tobs)
q0_hist <- matrix(0, imax, tmax)
q2_hist <- matrix(0, imax, tmax)
q4_hist <- matrix(0, imax, tmax)
f_hist <- matrix(0, imax, tmax)
n <- 1
for (k in 1:ntobs) {
  cat(sprintf("k= %d n = %d\n", k, n))
  if (n %in% tobs) {
    print("update ensemble")
    ens <- ensemble_update(xf, xe, hobs[, k], hmat, rmat)
    xf <- ens$xf
    xe <- ens$xe
  }
  if (k != ntobs) {
    print("run ensemble forecast")
    efcst <- run_ensemble(xf, xe, f, dtobs[k] + 1)
    xf <- efcst$xf
    xe <- efcst$xe
    f <- efcst$f
    q0_hist[, n:(dtobs[k] + n)] <- efcst$state$q0
    q2_hist[, n:(dtobs[k] + n)] <- efcst$state$q2
    q4_hist[, n:(dtobs[k] + n)] <- efcst$state$q4
    f_hist[, n:(dtobs[k] + n)] <- efcst$state$f
    n <- n + dtobs[k]
  } else {
    q0_hist[, n] <- xf[1:imax]
    q2_hist[, n] <- xf[1:imax + imax]
    q4_hist[, n] <- xf[1:imax + 2 * imax]
  }
}
state <- list(q0 = q0_hist, q2 = q2_hist, q4 = q4_hist, f = f_hist)
#saveRDS(state, "analysis_enkf.rds")
plot_waves(x, state)
