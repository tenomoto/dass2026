import numpy as np 
from numpy import cos, sin, sqrt
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

#tobs = np.array([0, 16, 32])
tobs = np.arange(0, tmax, 4)
ntobs = len(tobs)
s_obs = 0.1
htrue = eo.calc_h(true_state["q0"], true_state["q2"], true_state["q4"])
hobs = htrue[tobs, :] + rng.normal(0, s_obs, size = [ntobs, imax])

rmat = np.identity(imax) * s_obs ** 2

s_fgs = 1.0 # first guess standard deviation
s_sys = 0.0

q0 = sin(k0 * x)
q2 = sin(k2 * x)
q4 = sin(k4 * x)
xf = np.concatenate([q0, q2, q4])
pf = np.identity(xf.size) * s_fgs ** 2
f = true_state["f"]
qmat = np.identity(imax) * s_sys ** 2

guess_state = run.forward(q0, q2, q4, f, tmax)

def update(xf, pf, yo, rmat, infl = 0.1):
    pf *= (1.0 + infl)
    omax = len(yo)
    hxf = eo.calc_h(xf[0:imax], xf[imax:2 * imax], xf[2 * imax:3 * imax])
    hmat = np.zeros([omax, xf.size])
    hmat[:, 0:imax] = np.eye(omax)
    hmat[:, imax:2 * imax] = np.eye(omax)
    hmat[:, 2 * imax:3 * imax] = -1.0*np.eye(omax)
    d = yo - hxf
    kmat = pf @ hmat.T @ np.linalg.inv(hmat @ pf @ hmat.T + rmat)
    xf = xf + kmat @ d
    pf = (np.eye(pf.shape[0]) - kmat @ hmat) @ pf
    return {"xf": xf, "pf": pf}

dtobs = np.diff(tobs)
q0_hist = np.zeros([tmax, imax])
q2_hist = np.zeros([tmax, imax])
q4_hist = np.zeros([tmax, imax])
f_hist = np.zeros([tmax, imax])
q0_pf = np.zeros([tmax, imax])
q2_pf = np.zeros([tmax, imax])
q4_pf = np.zeros([tmax, imax])
n = 0
for k in range(ntobs):
    print(f"k= {k:d} n = {n:d}")
    if n == tobs[k]:
        print("update state")
        up = update(xf, pf, hobs[k, :], rmat)
        xf = up["xf"]
        pf = up["pf"]
        q0_hist[n, :] = xf[0:imax]
        q2_hist[n, :] = xf[imax:2 * imax]
        q4_hist[n, :] = xf[2 * imax:3 * imax]
        q0_pf[n, :] = pf[0:imax, imax//2]
        q2_pf[n, :] = pf[imax:2 * imax, imax//2]
        q4_pf[n, :] = pf[2 * imax:3 * imax, imax//2]
    if k < ntobs - 1:
        print("run forecast")
        nt = dtobs[k] + 1
        n_next = n + dtobs[k]
        f = true_state["f"][n:n_next + 1, :]
        fcst = run.forecast(xf, f, nt, prtb_f=True, qmat=qmat)
        xf = fcst["xf"]
        q0_hist[(n+1):(n_next + 1), :] = fcst["state"]["q0"][1:, :]
        q2_hist[(n+1):(n_next + 1), :] = fcst["state"]["q2"][1:, :]
        q4_hist[(n+1):(n_next + 1), :] = fcst["state"]["q4"][1:, :]
        f_hist[(n+1):(n_next + 1), :] = fcst["state"]["f"][1:, :]
        for i in range(1,nt):
            pf = run.lyapunov(pf, 1, qmat=qmat)
            q0_pf[n+i, :] = pf[0:imax, imax//2]
            q2_pf[n+i, :] = pf[imax:2 * imax, imax//2]
            q4_pf[n+i, :] = pf[2 * imax:3 * imax, imax//2]
        n = n_next

state = {"q0": q0_hist, "q2": q2_hist, "q4": q4_hist, "f": f_hist,
         "q0_pf": q0_pf, "q2_pf": q2_pf, "q4_pf": q4_pf}
run.print_mse(state, true_state)
plot_waves(x, state, true_state)
rmse = plot_rmse(state, true_state, guess_state)

vlim = np.max(np.abs(pf))
plt.matshow(pf, cmap='coolwarm', vmin=-vlim, vmax=vlim)
plt.colorbar()
plt.title('KF Pf')
plt.savefig('kf_pf.png', dpi=300)
plt.show()

fig, fig_s = plot_state(state)
fig.suptitle('KF')
fig.savefig('kf_state.png', dpi=300)
fig_s.suptitle('KF Pf')
fig_s.savefig('kf_state_pfrow.png', dpi=300)
plt.show()
