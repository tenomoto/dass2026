import numpy as np
import settings as st
import eqocean as eo

imax = st.imax
sgm = st.sgm
gma0 = st.gma0
tau = st.tau
rng = st.rng

tmax = 33

qi = rng.normal(size = imax)
df = np.zeros(imax) # tlm has not forcing
dh = np.zeros(imax) # no obs
eta = np.zeros(imax)
eps = 0

q = np.zeros([imax, tmax])
q[:, 0] = qi
for n in range(1, tmax):
    q[:, n] = eo.forward(q[:, n - 1], df, sgm, gma0, tau)

lhs = np.dot(q[:, -1], q[:, -1])
print(f"LHS: {lhs}")
p = q[:, -1]
for i in range(tmax - 1)[::-1]:
    p = eo.adjoint(p, dh, sgm, gma0, eps, tau)
rhs = np.dot(qi, p)
print(f"RHS: {rhs}")
print(f"LHS - RHS: {lhs - rhs:e}")

 
