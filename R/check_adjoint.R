source("settings.R")
source("eqocean.R")

seed <- 514
set.seed(seed)

tmax <- 33

qi <- rnorm(imax)
df <- numeric(imax) # tlm has not forcing
dh <- numeric(imax) # no obs
eta <- numeric(imax)
eps <- 0

q <- matrix(0, imax, tmax)
q[, 1] <- qi
for (n in 2:tmax) {
  q[, n] <- forward(q[, n - 1], df, sgm, gma0, tau)
}

lhs <- sum(q[, tmax] * q[, tmax])
print(paste("LHS:", lhs))
p <- q[, tmax]
for (i in (tmax - 1):1) {
  p <- adjoint(p, dh, sgm, gma0, eps, tau)
}
rhs <- sum(qi * p)
print(paste("RHS:", rhs))
print(paste("LHS - RHS:", lhs - rhs))

 
