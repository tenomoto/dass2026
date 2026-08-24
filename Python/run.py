import numpy as np
from numpy import sqrt

import settings as st
import eqocean as eo
import run

imax = st.imax
tmax = st.tmax
tau = st.tau
sgm = st.sgm
gma0, gma2, gma4, gmaf = st.gma0, st.gma2, st.gma4, st.gmaf
d0, d2, d4 = st.d0, st.d2, st.d4
eps0, eps2, eps4 = st.eps0, st.eps2, st.eps4

def forward(q0, q2, q4, f, nt, progf = True):
    q0_hist = np.zeros([nt, imax])
    q2_hist = np.zeros([nt, imax])
    q4_hist = np.zeros([nt, imax])
    q0_hist[0, :] = q0
    q2_hist[0, :] = q2
    q4_hist[0, :] = q4
    if progf: 
        for n in range(1, nt):
            f[n, :] = eo.forward(f[n - 1, :], 0, 0, gmaf, tau)
    for n in range(1, nt):
        q0 = eo.forward(q0, d0 * f[n, :], sgm, gma0, tau)
        q2 = eo.forward(q2, d2 * f[n, :], sgm, gma2, tau)
        q4 = eo.forward(q4, d4 * f[n, :], sgm, gma4, tau)
        q0_hist[n, :] = q0
        q2_hist[n, :] = q2
        q4_hist[n, :] = q4
    return {"q0":q0_hist, "q2":q2_hist, "q4":q4_hist, "f":f}

def tangent(q0, q2, q4, nt):
    q0_hist = np.zeros([nt, imax])
    q2_hist = np.zeros([nt, imax])
    q4_hist = np.zeros([nt, imax])
    q0_hist[0, :] = q0
    q2_hist[0, :] = q2
    q4_hist[0, :] = q4
    for n in range(1, nt):
        q0 = eo.tangent(q0, sgm, gma0, tau)
        q2 = eo.tangent(q2, sgm, gma2, tau)
        q4 = eo.tangent(q4, sgm, gma4, tau)
        q0_hist[n, :] = q0
        q2_hist[n, :] = q2
        q4_hist[n, :] = q4
    return {"q0":q0_hist, "q2":q2_hist, "q4":q4_hist}

def adjoint(dh, ds = None):
    p0 = np.zeros(imax) # tmax + 1
    p2 = np.zeros(imax)
    p4 = np.zeros(imax)
    if ds is None:
        ps = np.zeros([tmax, imax]) 
    for n in range(tmax)[::-1]:
        p0 = eo.adjoint(p0, dh[n, :], sgm, gma0, eps0, tau)
        p2 = eo.adjoint(p2, dh[n, :], sgm, gma2, eps2, tau)
        p4 = eo.adjoint(p4, dh[n, :], sgm, gma4, eps4, tau)
        if ds is not None and n > 0:
            ps[n, :] = -tau * (d0 * p0 + d2 * p2 + d4 * p4) - ds[n, :]
    if ds is None:
        return {"p0":p0, "p2":p2, "p4":p4}
    else:
        return {"p0":p0, "p2":p2, "p4":p4, "ps":ps}

def forecast(xf, f, nt, prtb_f=False, qmat=None):
    q0 = xf[0:imax]
    q2 = xf[imax:2 * imax]
    q4 = xf[2 * imax:3 * imax]
    if prtb_f and qmat is not None:
        f += st.rng.multivariate_normal(np.zeros(f.shape[1]), qmat, size=f.shape[0])
    state = run.forward(q0, q2, q4, f, nt)
    q0f = state["q0"][-1, :]
    q2f = state["q2"][-1, :]
    q4f = state["q4"][-1, :]
    xf = np.concatenate([q0f, q2f, q4f])
    return {"xf":xf, "state":state}

def ensemble(xf, xe, f, nt, prtb_f=False, qmat=None):
    nmem = xe.shape[1]
    q0m = xf[0:imax]
    q2m = xf[imax:(2 * imax)]
    q4m = xf[(2 * imax):(3 * imax)]
    q0_hist = np.zeros([nt, imax])
    q2_hist = np.zeros([nt, imax])
    q4_hist = np.zeros([nt, imax])
    f_hist = np.zeros([nt, imax])
    q0_sprd = np.zeros([nt, imax])
    q2_sprd = np.zeros([nt, imax])
    q4_sprd = np.zeros([nt, imax])
    f_sprd = np.zeros([nt, imax])
    for mem in range(nmem):
        q0 = q0m + xe[0:imax, mem]
        q2 = q2m + xe[imax:2 * imax, mem]
        q4 = q4m + xe[2 * imax:3 * imax, mem]
        fe = f.copy()
        if prtb_f and qmat is not None:
            fe += st.rng.multivariate_normal(np.zeros(f.shape[1]), qmat, size=f.shape[0])
        state = run.forward(q0, q2, q4, fe, nt)
        q0h, q2h, q4h, fh = state["q0"], state["q2"], state["q4"], state["f"]
        xe[:, mem] = np.concatenate([q0h[-1, :], q2h[-1, :], q4h[-1, :]])
        q0_hist += q0h
        q2_hist += q2h
        q4_hist += q4h
        f_hist += fh
        q0_sprd += q0h**2
        q2_sprd += q2h**2
        q4_sprd += q4h**2
        f_sprd += fh**2
    xm = np.mean(xe, axis = 1)
    xe = xe - xm[:, None]
    q0_hist = q0_hist / nmem
    q2_hist = q2_hist / nmem
    q4_hist = q4_hist / nmem
    f_hist = f_hist / nmem
    q0_sprd = q0_sprd / nmem - q0_hist**2
    q2_sprd = q2_sprd / nmem - q2_hist**2
    q4_sprd = q4_sprd / nmem - q4_hist**2
    f_sprd = f_sprd / nmem - f_hist**2
    q0_sprd = np.where(q0_sprd > 0, np.sqrt(q0_sprd), 0)
    q2_sprd = np.where(q2_sprd > 0, np.sqrt(q2_sprd), 0)
    q4_sprd = np.where(q4_sprd > 0, np.sqrt(q4_sprd), 0)
    f_sprd = np.where(f_sprd > 0, np.sqrt(f_sprd), 0)
    mstate = {"q0":q0_hist, "q2":q2_hist, "q4":q4_hist, "f":f_hist,
              "q0s":q0_sprd, "q2s":q2_sprd, "q4s":q4_sprd, "fs":f_sprd}
    f = mstate["f"]
    return {"xf":xm, "xe":xe, "f":f, "state":mstate}

# Pf = M @ Pa @ M.T + Q
def lyapunov(pa, nt, qmat=None):
    pf = pa.copy()
    for i in range(pf.shape[0]):
        e = np.zeros(pf.shape[0])
        e[i] = 1.0
        for n in range(nt)[::-1]:
            p0 = eo.adjoint(e[:imax], 0.0, sgm, gma0, eps0, tau)
            p2 = eo.adjoint(e[imax:2 * imax], 0.0, sgm, gma2, eps2, tau)
            p4 = eo.adjoint(e[2 * imax:3 * imax], 0.0, sgm, gma4, eps4, tau)
            e = np.concatenate([p0, p2, p4])
        e = pa @ e
        for n in range(nt):
            p0 = eo.tangent(e[:imax], sgm, gma0, tau)
            p2 = eo.tangent(e[imax:2 * imax], sgm, gma2, tau)
            p4 = eo.tangent(e[2 * imax:3 * imax], sgm, gma4, tau)
            e = np.concatenate([p0, p2, p4])
        pf[:, i] = e
    if qmat is not None:
        gmat = np.zeros([pf.shape[0], qmat.shape[0]])
        gmat[0:imax, :] = d0*np.eye(imax)
        gmat[imax:2 * imax, :] = d2*np.eye(imax)
        gmat[2 * imax:3 * imax, :] = d4*np.eye(imax)
        pf += gmat @ qmat @ gmat.T
    return pf

def calc_cost(dh, ds, sh = 1, ss = 1):
    return 0.5 * (np.sum(dh * dh) / sh ** 2 + np.sum(ds * ds) / ss ** 2)

def print_diag(epoch, cost, gnorm, message = ""): 
    print(f"epoch: {epoch:d}, cost={cost:f} gnorm={gnorm:e} {message}")

def print_mse(state, tstate):
    ha = eo.calc_h(state["q0"][-1, :], state["q2"][-1, :], state["q4"][-1, :])
    ht = eo.calc_h(tstate["q0"][-1, :], tstate["q2"][-1, :], tstate["q4"][-1, :])
    e = sum((ha - ht) ** 2)
    print(f"MSE = {e:e}")
