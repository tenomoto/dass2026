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

tobs <- c(0, 16, 32) + 1
ntobs <- length(tobs)
s_obs <- 0
htrue <- calc_h(true_state$q0, true_state$q2, true_state$q4)
hobs <- htrue[, tobs] + matrix(rnorm(imax * ntobs, 0, s_obs), imax, ntobs)

rmat <- diag(imax) * s_obs^2

nmem <- 200

s_ens <- 0.01
q0 <- sin(k0 * x)
q2 <- sin(k2 * x)
q4 <- sin(k4 * x)
xf <- c(q0, q2, q4)
xe <- matrix(rnorm(3 * imax * nmem, 0, s_ens), 3 * imax, nmem)
f <- true_state$f

guess_state <- run_forward(q0, q2, q4, f, tmax)

ensemble_update <- function(xf, xe, yo, rmat, infl = 0.1, std = 0.01) {
  nmem <- ncol(xe)
  omax <- length(yo)
  xe <- sqrt(1 + infl) * xe
  hxe <- matrix(0, omax, nmem)
  for (mem in 1:nmem) {
    q0 <- xf[1:imax] + xe[1:imax, mem]
    q2 <- xf[1:imax + imax] + xe[1:imax + imax, mem]
    q4 <- xf[1:imax + 2 * imax] + xe[1:imax + 2 * imax, mem]
    hxe[, mem] <- calc_h(q0, q2, q4)
  }
  hxm <- rowMeans(hxe)
  hxe <- hxe - hxm
  # K = X^f*Y^T*(Y*Y^T + (N-1)*R)^{-1}
  kmat <- xe %*% t(hxe) %*% solve(hxe %*% t(hxe) + (nmem - 1) * rmat)
  xf <- xf + kmat %*% (yo - hxm)
  ye <- matrix(rnorm(omax * nmem, 0, std), omax, nmem)
  xe <- xe + kmat %*% (ye - hxe)
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
  if (n == tobs[k]) {
    print("update ensemble")
    ens <- ensemble_update(xf, xe, hobs[, k], rmat)
    xf <- ens$xf
    xe <- ens$xe
    q0_hist[, n] <- xf[1:imax]
    q2_hist[, n] <- xf[1:imax + imax]
    q4_hist[, n] <- xf[1:imax + 2 * imax]
  }
  if (k < ntobs) {
    print("run ensemble forecast")
    n_next <- n + dtobs[k]
    nt <- dtobs[k] + 1
    f <- true_state$f[, n:n_next]
    efcst <- run_ensemble(xf, xe, f, dtobs[k] + 1)
    xf <- efcst$xf
    xe <- efcst$xe
    q0_hist[, (n + 1):n_next] <- efcst$state$q0[, 2:nt]
    q2_hist[, (n + 1):n_next] <- efcst$state$q2[, 2:nt]
    q4_hist[, (n + 1):n_next] <- efcst$state$q4[, 2:nt]
    f_hist[, (n + 1):n_next] <- efcst$state$f[, 2:nt]
    n <- n_next
  }
}
state <- list(q0 = q0_hist, q2 = q2_hist, q4 = q4_hist, f = f_hist)
print_mse(state, true_state)
plot_waves(x, state, true_state)
