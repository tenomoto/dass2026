source("settings.R")
source("eqocean.R")

run_forward <- function(q0, q2, q4, fin, tmax) {
  q0_hist <- matrix(0, imax, tmax)
  q2_hist <- matrix(0, imax, tmax)
  q4_hist <- matrix(0, imax, tmax)
  q0_hist[, 1] <- q0
  q2_hist[, 1] <- q2
  q4_hist[, 1] <- q4
  if (is.vector(fin)) {
    f_hist <- matrix(0, imax, tmax)
    f_hist[, 1] <- fin
    f <- fin
  }
  for (n in 2:tmax) {
    if (is.vector(fin)) {
      f <- forward(f, 0, 0, gmaf, tau)
    } else {
      f <- fin[, n]
    }
    q0 <- forward(q0, d0 * f, sgm, gma0, tau)
    q2 <- forward(q2, d2 * f, sgm, gma2, tau)
    q4 <- forward(q4, d4 * f, sgm, gma4, tau)
    q0_hist[, n] <- q0
    q2_hist[, n] <- q2
    q4_hist[, n] <- q4
    if (is.vector(fin)) f_hist[, n] <- f 
  }
  if (is.matrix(fin)) {
    f_hist <- fin
  }
  list(q0 = q0_hist, q2 = q2_hist, q4 = q4_hist, f = f_hist)
}

run_adjoint <- function(dh, ds = NULL) {
  p0 <- numeric(imax) # tmax + 1
  p2 <- numeric(imax)
  p4 <- numeric(imax)
  if (!is.null(ds)) ps <- matrix(0, imax, tmax)
  for (n in tmax:1) {
    p0 <- adjoint(p0, dh[, n], sgm, gma0, eps0, tau)
    p2 <- adjoint(p2, dh[, n], sgm, gma2, eps2, tau)
    p4 <- adjoint(p4, dh[, n], sgm, gma4, eps4, tau)
    if (!is.null(ds) & n > 1)  ps[, n] <- -tau * (d0 * p0 + d2 * p2 + d4 * p4) - ds[, n]
  }
  if (is.null(ds)) {
    list(p0 = p0, p2 = p2, p4 = p4)
  } else {
    list(p0 = p0, p2 = p2, p4 = p4, ps = ps)
  }
}

run_ensemble <- function(xf, xe, f, nt) {
  # input:
  #  xf: ensemble mean
  #  xe: ensemble perturbation
  #  f : forcing
  #  nt: number of steps
  q0m <- xf[1:imax]
  q2m <- xf[1:imax + imax]
  q4m <- xf[1:imax + 2 * imax]
  q0_hist <- matrix(0, imax, nt)
  q2_hist <- matrix(0, imax, nt)
  q4_hist <- matrix(0, imax, nt)
  f_hist <- matrix(0, imax, nt)
  q0s_hist <- matrix(0, imax, nt)
  q2s_hist <- matrix(0, imax, nt)
  q4s_hist <- matrix(0, imax, nt)
  for (mem in 1:nmem) {
    q0 <- q0m + xe[1:imax, mem]
    q2 <- q2m + xe[1:imax + imax, mem]
    q4 <- q4m + xe[1:imax + 2 * imax, mem]
    state <- run_forward(q0, q2, q4, f, nt)
    xe[, mem] <- c(state$q0[, nt], state$q2[, nt], state$q4[, nt])
    q0_hist <- q0_hist + state$q0
    q2_hist <- q2_hist + state$q2
    q4_hist <- q4_hist + state$q4
    f_hist <- f_hist + state$f
    q0s_hist <- q0s_hist + state$q0^2
    q2s_hist <- q2s_hist + state$q2^2
    q4s_hist <- q4s_hist + state$q4^2
  }
  xm <- rowMeans(xe)
  xe <- xe - xm
  q0_hist <- q0_hist / nmem
  q2_hist <- q2_hist / nmem
  q4_hist <- q4_hist / nmem
  q0s_hist <- (q0s_hist/nmem) - q0_hist^2
  q2s_hist <- (q2s_hist/nmem) - q2_hist^2
  q4s_hist <- (q4s_hist/nmem) - q4_hist^2
  mstate <- list(q0 = q0_hist, q2 = q2_hist, q4 = q4_hist,
                 f = f_hist / nmem,
                 q0s = sqrt(q0s_hist), q2s = sqrt(q2s_hist), q4s = sqrt(q4s_hist))
  f <- mstate$f[, nt]
  list(xf = xm, xe = xe, f = f, state = mstate)
}

calc_cost <- function(dh, ds, sh = 1, ss = 1) 0.5 * (sum(dh^2) / sh^2 + sum(ds^2) / ss^2)

print_diag <- function(epoch, cost, gnorm, message = "") {
  cat(sprintf("epoch: %d, cost=%f gnorm=%e %s\n", epoch, cost, gnorm, message))
}

print_mse <- function(state, tstate) {
  ha <- calc_h(state$q0[,tmax], state$q2[, tmax], state$q4[, tmax])
  ht <- calc_h(tstate$q0[,tmax], tstate$q2[, tmax], tstate$q4[, tmax])
  cat(sprintf("MSE = %e\n", sum((ha - ht)^2)))
}
