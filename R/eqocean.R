cshift <- function(n, k) {
  ((1:n + k - 1) %% n + 1)
}

forward <- function(q, df, sgm, gma, tau) {
  n <- length(q)
  (q - 0.5 * gma * (q[cshift(n, 1)] - q[cshift(n, -1)]) +
    　 0.5 * gma^2 *(q[cshift(n, 1)] - 2 * q + q[cshift(n, -1)]) +
       tau * df) / (1 + sgm * tau)
}

adjoint <- function(p, q0, q2, q4, eta, sgm, gma, epsh) {
  (p - 0.5 * gma * (p[cshift(n, -1)] - p[cshift(n, 1)]) +
       0.5 * gma^2 (p[cshift(n, -1)] - 2 * p + p[cshift(n, 1)]) +
       epsh * (eta - q0 - q2 + q4)
  ) / (1 + sgm * tau)
}

retrieve <- function(tdp, phi, fflg) {
  (-tdp - phi) / fflg
}
