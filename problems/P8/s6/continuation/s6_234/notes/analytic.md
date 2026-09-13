# First sheet and the complete complex-disk bound

The mixed bubble uses the full Feynman denominator

Q(x,s)=xM_H²+(1-x)-x(1-x)s.

For nonreal s and0<x<1 its imaginary part is nonzero. At the physical threshold,

Q(x,(M_H+1)²)=[(M_H+1)x-1]².

For smaller real s, Q is positive on the full parameter interval. At the pseudothreshold it is

Q(x,(M_H-1)²)=[(M_H-1)x+1]²>0.

Thus the mixed first bubble's first sheet is cut along[(M_H+1)²,infinity); the pseudothreshold is not an extra first-sheet cut. Higher-loop cuts containing only light particles are a different obligation and can begin below this threshold.

For the disk centered at s1, write d=s-1 and F=Q(x,1). The exact positive identity

(1-x)/M_H²-alpha=(1-x)³/(M_H² F)

implies0<=alpha<=(1-x)/M_H². If abs(d)<=R and r=R/M_H²<1, then abs(alpha d)<=r. The analytic logarithm remainder has the convergent series

-Log(1-u)-u=sum_(n>=2) u^n/n,

and therefore its modulus is at most abs(u)²/[2(1-r)]. This follows by absolute convergence and1/n<=1/2; it applies to complex u, not just a sampled positive real ray. Integrating the squared alpha majorant gives

integral alpha²<=1/(3M_H4),
abs(Pi_OS)<=g² abs(d)²/[96pi²M_H4(1-r)].

These dominations also justify differentiation under the integral. The quotient q=Pi_OS/d extends holomorphically across d0, with q(0)=0.

At R=10^196, r<r0=32/625. Using pi²>9,

abs(q)<=epsilon,
epsilon=g² r0/[864 M_H²(1-r0)].

Exact arithmetic gives0<epsilon<10^-209 and epsilon/(1-epsilon)<10^-209. In particular1+q never vanishes on the entire disk. The reciprocal of the ONE-LOOP-TRUNCATED denominator is

1/[d+Pi_OS]=1/[d(1+q)].

It has only the anchor simple pole there, with unit residue. The continuously extended propagator ratio differs from1 by at most epsilon/(1-epsilon). Both conclusions follow on the complete complex disk; a finite contour sample is not the proof.

This reciprocal is an algebraic diagnostic of the truncated inverse. Its implicit higher powers of Pi_OS do not constitute a complete higher-loop resummation and do not bound the true propagator. It gives no full angular four-point quantum error, physical cutoff, exact pole exclusion outside the disk, bounce stability or original affine result.
