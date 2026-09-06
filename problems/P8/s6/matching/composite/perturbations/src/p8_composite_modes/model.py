"""Positive-root FLRW symbols, in the frozen +--- physical composite frame."""

import sympy as sp

a, Ng, y, c, G, F, M4 = sp.symbols("a Ng y c G F m4", positive=True)
alpha, beta = sp.symbols("alpha beta", positive=True)
rho, p = sp.symbols("rho p", real=True)
betas = sp.symbols("beta0:5", real=True)
r = alpha+beta*y
s = alpha+beta*c
P = M4*(betas[1]+2*betas[2]*y+betas[3]*y**2)
n = rho+p
b = a*y
Nf = Ng*c
Ae = a*r
Ne = Ng*s


def cancel(value):
    return sp.factor(sp.cancel(value))
