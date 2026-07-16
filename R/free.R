source("settings.R")
source("eqocean.R")

tmax <- 32

q0 <- cos(k0 * x)
q2 <- cos(k2 * x)
q4 <- cos(k4 * x)
f0 <- 0.1 * cos(kf * x)
f <- 0.1 * cos(kf * x)

q0_hist <- matrix(0, imax, tmax)
q2_hist <- matrix(0, imax, tmax)
q4_hist <- matrix(0, imax, tmax)
f_hist <- matrix(0, imax, tmax)
q0_hist[, 1] <- q0
q2_hist[, 1] <- q2
q4_hist[, 1] <- q4
f_hist[, 1] <- f 
for (i in 2:tmax) {
  q0 <- forward(q0, d0 * f, sgm, gma0, tau)
  q2 <- forward(q2, d2 * f, sgm, gma2, tau)
  q4 <- forward(q4, d4 * f, sgm, gma4, tau)
  f <- forward(f, 0, 0, gmaf, tau)
  plot_waves(x, q0, q2, q4, f, i)
  q0_hist[, i] <- q0
  q2_hist[, i] <- q2
  q4_hist[, i] <- q4
  f_hist[, i] <- f 
}
hist <- list(q0 = q0_hist, q2 = q2_hist, q4 = q4_hist, f = f_hist)
saveRDS(hist, file = "free.rds")
