cshift <- function(n, k) {
  ((1:n + k - 1) %% n + 1)
}

forward <- function(q, df, sgm, gma, tau) {
  n <- length(q)
  (q - 0.5 * gma * (q[cshift(n, 1)] - q[cshift(n, -1)]) +
    　 0.5 * gma^2 *(q[cshift(n, 1)] - 2 * q + q[cshift(n, -1)]) +
       tau * df) / (1 + sgm * tau)
}

adjoint <- function(p, h, eta, sgm, gma, epsh, tau) {
  n <- length(p)
  (p - 0.5 * gma * (p[cshift(n, -1)] - p[cshift(n, 1)]) +
       0.5 * gma^2 * (p[cshift(n, -1)] - 2 * p + p[cshift(n, 1)]) +
       epsh * (eta - h)
  ) / (1 + sgm * tau)
}

calc_h <- function(q0, q2, q4) q0 + q2 - q4

retrieve <- function(tdp, phi, fflg) {
  ifelse(fflg == 1, -tdp - phi, 0)
}
