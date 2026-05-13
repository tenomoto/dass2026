source("eqocean.R")

d0 <- 1
d2 <- 0.5
d4 <- 0.125

tau <- 1
sgm <- 0.05
tmax <- 32
dx <- 1
imax <- 75
x <- 2 * pi / imax * (1:imax - 1)
k0 <- 1
k2 <- 3
k4 <- 7
kf <- 2
c0 <- 1
c2 <- -c0 / 3
c4 <- -c0 / 7
cf <- 0.5 * c0
gma0 <- c0 * tau / dx
gma2 <- c2 * tau / dx
gma4 <- c4 * tau / dx
gmaf <- cf * tau / dx

q0 <- cos(k0 * x)
q2 <- cos(k2 * x)
q4 <- cos(k4 * x)
f0 <- 0.1 * cos(kf * x)
f <- 0.1 * cos(kf * x)

for (i in 1:(tmax - 1)) {
  q0 <- forward(q0, d0 * f, sgm, gma0, tau)
  q2 <- forward(q2, d2 * f, sgm, gma2, tau)
  q4 <- forward(q4, d4 * f, sgm, gma4, tau)
  f <- forward(f, 0, 0, gmaf, tau)
}
h <- q0 + q2 - q4
qmax <- max(c(q0, q2, q4, h, f))
qmin <- min(c(q0, q2, q4, h, f))
plot(x, q0, type = "l", ylim = c(qmin, qmax), col = "blue",
     main = "equatorial ocean waves",
     xlab = "phase", ylab = "amplitude", xaxt = "n")
axis(1, c(0, pi/2, pi, 3 * pi / 2, 2 * pi),
     c("0", "π/2", "π", "3π/2", "2π"))
abline(h = 0, lty = 2)
lines(x, q2, col = "red")
lines(x, q4, col = "orange")
lines(x, h, col = "black")
lines(x, f, col = "purple")
legend("topright", c("q0", "q2", "q4", "h", "f"), lty = 1,
       col = c("blue", "red", "orange", "black", "purple"))

