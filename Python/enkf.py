import numpy as np
from numpy import cos, sin, sqrt

import settings as st
import eqocean as eo
import run
from plot import plot_waves, plot_rmse

x = st.x
tmax = st.tmax
imax = st.imax
k0, k2, k4, kf = st.k0, st.k2, st.k4, st.kf
rng = st.rng

q0 = cos(k0 * x)
q2 = cos(k2 * x)
q4 = cos(k4 * x)
f = np.zeros([tmax, imax])
f[0, :] = 0.1 * cos(kf * x)

true_state = run.forward(q0, q2, q4, f, tmax)

tobs = np.array([0, 16, 32])
ntobs = len(tobs)
s_obs = 0.01
htrue = eo.calc_h(true_state["q0"], true_state["q2"], true_state["q4"])
hobs = htrue[tobs, :] + rng.normal(0, s_obs, size = [ntobs, imax])

rmat = np.identity(imax) * s_obs ** 2

nmem = 100

s_mod = 0.01
s_ens = 1.0

q0 = sin(k0 * x)
q2 = sin(k2 * x)
q4 = sin(k4 * x)
xf = np.concatenate([q0, q2, q4])
xe = rng.normal(0, s_ens, size = [3 * imax, nmem])
f = true_state["f"]

guess_state = run.forward(q0, q2, q4, f, tmax)

def ensemble_update(xf, xe, yo, rmat, infl = 0.1, std = 0.01):
    nmem = xe.shape[1]
    omax = len(yo)
    hxe = np.zeros([omax, nmem])
    for mem in range(nmem):
        q0 = xf[0:imax] + xe[0:imax, mem]
        q2 = xf[imax:2 * imax] + xe[imax:2 * imax, mem]
        q4 = xf[2 * imax:3 * imax] + xe[2 * imax:3 * imax, mem]
        hxe[:, mem] = eo.calc_h(q0, q2, q4)
    hxm = np.mean(hxe, axis = 1)
    hxe = hxe - hxm[:, None]
    xe = sqrt(1 + infl) * xe
    kmat = xe @ hxe.T @ np.linalg.inv(hxe @ hxe.T + (nmem - 1) * rmat)
    xf = xf + kmat @ (yo - hxm)
    ye = rng.normal(0, std, size = [omax, nmem])
    xe = xe + kmat @ (ye - hxe)
    return {"xf": xf, "xe": xe}

dtobs = np.diff(tobs)
q0_hist = np.zeros([tmax, imax])
q2_hist = np.zeros([tmax, imax])
q4_hist = np.zeros([tmax, imax])
f_hist = np.zeros([tmax, imax])
n = 0
for k in range(ntobs):
    print(f"k= {k:d} n = {n:d}")
    if n in tobs:
        print("update ensemble")
        ens = ensemble_update(xf, xe, hobs[k, :], rmat)
        xf = ens["xf"]
        xe = ens["xe"]
    if k < ntobs - 1:
        print("run ensemble forecast")
        f = true_state["f"][n:(dtobs[k] + n + 1), :]
        efcst = run.ensemble(xf, xe, f, dtobs[k] + 1)
        xf = efcst["xf"]
        xe = efcst["xe"]
        f = efcst["f"]
        q0_hist[n:(dtobs[k] + n + 1), :] = efcst["state"]["q0"]
        q2_hist[n:(dtobs[k] + n + 1), :] = efcst["state"]["q2"]
        q4_hist[n:(dtobs[k] + n + 1), :] = efcst["state"]["q4"]
        f_hist[n:(dtobs[k] + n + 1), :] = efcst["state"]["f"]
        n = n + dtobs[k]
    else:
        q0_hist[n, :] = xf[0:imax]
        q2_hist[n, :] = xf[imax:2 * imax]
        q4_hist[n, :] = xf[2 * imax:3 * imax]

state = {"q0": q0_hist, "q2": q2_hist, "q4": q4_hist, "f": f_hist}
run.print_mse(state, true_state)
plot_waves(x, state, true_state)
rmse = plot_rmse(state, true_state, guess_state)
