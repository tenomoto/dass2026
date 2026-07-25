import numpy as np
from numpy import cos, sin, fabs
import settings as st
import eqocean as eo
import run
from plot import plot_waves

x = st.x
tmax = st.tmax
imax = st.imax
k0, k2, k4, kf = st.k0, st.k2, st.k4, st.kf
rng = st.rng

q0 = cos(k0 * x)
q2 = cos(k2 * x)
q4 = cos(k4 * x)
f = 0.1 * cos(kf * x)

true_state = run.forward(q0, q2, q4, f, tmax)

tobs = np.array([0, 16, 32])
ntobs = len(tobs)
sobs = 0
htrue = eo.calc_h(true_state["q0"], true_state["q2"], true_state["q4"])
hobs = htrue[tobs, :] + rng.normal(0, sobs, [ntobs, imax])

niter = 1000
ftol = 1e-15
ctol = 1e-15
gtol = 1e-15

smod = 0.01
q0 = sin(k0 * x)
q2 = sin(k2 * x)
q4 = sin(k4 * x)
f = 0.1 * cos(kf * x)

#guess_state = run.forward(q0, q2, q4, f, tmax)

lr = 0.1
ocost = 1e8
ognorm = 1
fcg = 1 # set to 1 to use CG
d = np.zeros(3 * imax)

dh = np.zeros([tmax, imax])
state = run.forward(q0, q2, q4, f, tmax)
for ep in range(niter):
    state = run.forward(q0, q2, q4, f, tmax)
    h = eo.calc_h(state["q0"], state["q2"], state["q4"])
    dh[:, :] = 0
    dh[tobs, :] = hobs - h[tobs, :]
    adj = run.adjoint(dh)
    cost = run.calc_cost(dh, 0)
    xf = np.concatenate([q0, q2, q4])
    g = -np.concatenate([adj["p0"], adj["p2"], adj["p4"]])
    gnorm = np.dot(g, g)
    bta = gnorm / ognorm
    d = -g + fcg * bta * d
    xa = xf + lr * d
    q0 = xa[0:imax]
    q2 = xa[imax:(2 * imax)]
    q4 = xa[(2 * imax):(3 * imax)]
    dcost = fabs(cost - ocost)
    if ep % 10 == 0:
        run.print_diag(ep, cost, gnorm)
    if cost < ftol:
        run.print_diag(ep, cost, gnorm, "stopping due to ftol")
        break
    if dcost < ctol: 
        run.print_diag(ep, cost, gnorm, "stopping due to ctol")
        break
    if gnorm < gtol:
        run.print_diag(ep, cost, gnorm, "stopping due to gtol")
        break
    ocost = cost
    ognorm = gnorm
state = run.forward(q0, q2, q4, f, tmax)
run.print_mse(state, true_state)
#saveRDS(state, "analysis_var.rds")
#plot_rmse(state, true_state, guess_state)
plot_waves(x, state, true_state)
#plot_diff(x, state, true_state)
