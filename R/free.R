source("settings.R")
source("eqocean.R")
source("plot.R")

tmax <- 32

q0 <- cos(k0 * x)
q2 <- cos(k2 * x)
q4 <- cos(k4 * x)
f0 <- 0.1 * cos(kf * x)
f <- 0.1 * cos(kf * x)

plot_waves(x, q0, q2, q4, f, i)
for (i in 2:tmax) {
  q0 <- forward(q0, d0 * f, sgm, gma0, tau)
  q2 <- forward(q2, d2 * f, sgm, gma2, tau)
  q4 <- forward(q4, d4 * f, sgm, gma4, tau)
  f <- forward(f, 0, 0, gmaf, tau)
  plot_waves(x, q0, q2, q4, f, i)
}
