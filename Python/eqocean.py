import numpy as np

def cshift(x, k):
    return np.roll(x, -k)

def forward(q, df, sgm, gma, tau): 
    return (q - 0.5 * gma * (cshift(q, 1) - cshift(q, -1)) +
            0.5 * gma ** 2 *(cshift(q, 1) - 2 * q + cshift(q, -1)) +
            tau * df) / (1 + sgm * tau)

def tangent(q, sgm, gma, tau):
    return (q - 0.5 * gma * (cshift(q, 1) - cshift(q, -1)) + 
            0.5 * gma ** 2 *(cshift(q, 1) - 2 * q + cshift(q, -1))) / (1 + sgm * tau)

def adjoint(p, dh, sgm, gma, eps, tau):
    return (p - 0.5 * gma * (cshift(p, -1) - cshift(p, 1)) +
            0.5 * gma ** 2 * (cshift(p, -1) - 2 * p + cshift(p, 1)) +
            eps * dh) / (1 + sgm * tau)

def calc_h(q0, q2, q4):
    return q0 + q2 - q4
