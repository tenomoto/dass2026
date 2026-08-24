import numpy as np
from numpy import cos, sin, sqrt, pi
import matplotlib.pyplot as plt
import settings as st
import eqocean as eo
import run
from plot import plot_waves, plot_rmse, plot_state

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
fig = plot_state(true_state)
fig.suptitle('True state')
fig.savefig('true_state.png', dpi=300)
plt.show(block=False)
plt.pause(1.0)
plt.close()

#tobs = np.array([0, 16, 32])
tobs = np.arange(0, tmax, 4)
ntobs = len(tobs)
s_obs = 0.1
htrue = eo.calc_h(true_state["q0"], true_state["q2"], true_state["q4"])
hobs = htrue[tobs, :] + rng.normal(0, s_obs, size = [ntobs, imax])
for to in tobs:
    plt.figure(figsize=(8,4),constrained_layout=True)
    plt.plot(x, htrue[to, :], label = 'truth')
    plt.plot(x, hobs[np.where(tobs == to)[0][0], :], 'x', lw=0.0, label = 'obs')
    plt.xlabel('Location')
    plt.ylabel('h')
    plt.title(f'Observations at t = {to:d}')
    plt.legend()
    plt.savefig(f'obs_t{to:d}.png', dpi=300)
    plt.show(block=False)
    plt.pause(1.0)
    plt.close()

rmat = np.identity(imax) * s_obs ** 2

nmem = 200

s_ens = 0.4
s_sys = 0.0

q0 = sin(k0 * x)
q2 = sin(k2 * x)
q4 = sin(k4 * x)
xf = np.concatenate([q0, q2, q4])
xe = rng.normal(0, s_ens, size = [3 * imax, nmem])
#q0e = sin(k0 * x[:, None] + rng.normal(0, 2.0*pi, size = [imax, nmem]))
#q2e = sin(k2 * x[:, None] + rng.normal(0, 2.0*pi, size = [imax, nmem]))
#q4e = sin(k4 * x[:, None] + rng.normal(0, 2.0*pi, size = [imax, nmem]))
#xe = np.concatenate([q0e, q2e, q4e], axis=0)*s_ens
#xe -= np.mean(xe, axis=1)[:, None]
f = true_state["f"]
qmat = np.identity(imax) * s_sys ** 2

guess_state = run.forward(q0, q2, q4, f, tmax)
fig = plot_state(guess_state)
fig.suptitle('First guess')
fig.savefig('guess_state.png', dpi=300)
plt.show(block=False)
plt.pause(1.0)
plt.close()

def ensemble_update(xf, xe, yo, rmat, infl = 0.1, std = s_obs):
    nmem = xe.shape[1]
    omax = len(yo)
    xe = sqrt(1 + infl) * xe
    hxe = np.zeros([omax, nmem])
    for mem in range(nmem):
        q0 = xf[0:imax] + xe[0:imax, mem]
        q2 = xf[imax:2 * imax] + xe[imax:2 * imax, mem]
        q4 = xf[2 * imax:3 * imax] + xe[2 * imax:3 * imax, mem]
        hxe[:, mem] = eo.calc_h(q0, q2, q4)
    hxm = np.mean(hxe, axis = 1)
    hxe = hxe - hxm[:, None]
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
q0_sprd = np.zeros([tmax, imax])
q2_sprd = np.zeros([tmax, imax])
q4_sprd = np.zeros([tmax, imax])
f_sprd = np.zeros([tmax, imax])
n = 0
for k in range(ntobs):
    print(f"k= {k:d} n = {n:d}")
    if n == tobs[k]:
        print("update ensemble")
        ens = ensemble_update(xf, xe, hobs[k, :], rmat)
        xf = ens["xf"]
        xe = ens["xe"]
        xall = xf[:, None] + xe
        q0_hist[n, :] = xf[0:imax]
        q2_hist[n, :] = xf[imax:2 * imax]
        q4_hist[n, :] = xf[2 * imax:3 * imax]
        q0_sprd[n, :] = np.std(xall[0:imax, :], axis=1)
        q2_sprd[n, :] = np.std(xall[imax:2 * imax, :], axis=1)
        q4_sprd[n, :] = np.std(xall[2 * imax:3 * imax, :], axis=1)
    if k < ntobs - 1:
        print("run ensemble forecast")
        nt = dtobs[k] + 1
        n_next = n + dtobs[k]
        f = true_state["f"][n:n_next + 1, :]
        efcst = run.ensemble(xf, xe, f, nt, prtb_f=True, qmat=qmat)
        xf = efcst["xf"]
        xe = efcst["xe"]
        q0_hist[(n + 1):(n_next + 1), :] = efcst["state"]["q0"][1:, :]
        q2_hist[(n + 1):(n_next + 1), :] = efcst["state"]["q2"][1:, :]
        q4_hist[(n + 1):(n_next + 1), :] = efcst["state"]["q4"][1:, :]
        f_hist[(n + 1):(n_next + 1), :] = efcst["state"]["f"][1:, :]
        q0_sprd[(n + 1):(n_next + 1), :] = efcst["state"]["q0s"][1:, :]
        q2_sprd[(n + 1):(n_next + 1), :] = efcst["state"]["q2s"][1:, :]
        q4_sprd[(n + 1):(n_next + 1), :] = efcst["state"]["q4s"][1:, :]
        f_sprd[(n + 1):(n_next + 1), :] = efcst["state"]["fs"][1:, :]
        n = n_next

state = {"q0": q0_hist, "q2": q2_hist, "q4": q4_hist, "f": f_hist, "q0s": q0_sprd, "q2s": q2_sprd, "q4s": q4_sprd, "fs": f_sprd}
run.print_mse(state, true_state)
plot_waves(x, state, true_state)
rmse = plot_rmse(state, true_state, guess_state)
fig, fig_s = plot_state(state)
fig.suptitle('EnKF')
fig.savefig('enkf_state.png', dpi=300)
fig_s.suptitle('EnKF spread')
fig_s.savefig('enkf_state_sprd.png', dpi=300)
plt.show()