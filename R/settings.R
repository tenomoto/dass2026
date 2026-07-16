d0 <- 1
d2 <- 0.5
d4 <- 0.125

tau <- 1
sgm <- 0.05
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
eps0 <- 1
eps2 <- 1
eps4 <- -1
