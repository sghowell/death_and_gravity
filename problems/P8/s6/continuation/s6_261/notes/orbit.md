# Complete density-weight gauge orbit

Write gamma=a^2 exp(2v) exp(h), Tr h=0, and
Q=(det gamma)^(1/3) gamma^-1=exp(-h). Its S258 gauge is chi=div Q=0.
For w near 2/3 define Q_w=(det gamma)^(w/2) gamma^-1. With d=w-2/3,

Q_w=exp[3d(v+log a)] Q.

Differentiate BEFORE restriction:
partial_w chi_w=3(v+log a)chi+3Q grad v.
The first term is present off the slice. On chi=0 the infinitesimal
coordinate reconstruction solves

M_Q xi=-3Q grad v,
M_Q xi=-Q^(jk)partial_j partial_k xi-(Q grad div xi)/3.

This is the full gauge-surface operator derived from the density Lie
variation; its off-slice terms remain in S258/S260. The local mean-zero
complement and translations are treated exactly as there. Since div Q=0,
Q grad v=div(vQ) has zero periodic integral. Coercivity and elliptic
regularity apply in the already stated local shape domain, not on a
general global metric space and not with an infinite-volume uniform gap.

The logarithmic volume and scalar lapse vary as
delta v=xi dot grad v+div xi/3 and delta N=xi dot grad N.
The first-order nonzero-mode equation at Q=I gives
xi1(k)=-(9/4)i k v(k)/k^2 and delta v1=3v_perp/4.
The zero mode is not divided by k^2.

## An exact finite periodic family

Take a 2pi-periodic coordinate y and retain two arbitrary smooth tensor
profiles t(y),z(y). The transverse matrix

B=[[exp(t),z],[z,(1+z^2)exp(-t)]]

is positive definite and has determinant one. Use
gamma=exp(2v(y))diag(1,B(y)). This includes two tensor profiles and is
exactly on the original three-component divergence slice, not merely
linear TT. For w<2 set

beta=(3w-2)/(2-w), K=<exp(-beta v)>_y,
dx/dy=exp(-beta v)/K, y=phi(x), phi'=K exp(beta v(phi)).

The strictly positive normalized derivative defines an orientation-
preserving circle diffeomorphism. A residual constant translation can
fix the origin or displacement convention; spatial averages are
unaffected. The full pulled metric is
exp(2v(phi))diag(phi'^2,B(phi)). Direct inversion gives
Q_w^(xx)=K^(w-2) and Q_w^(xy)=Q_w^(xz)=0. Hence all three gauge
divergences vanish exactly; no tensor profile has been dropped.
This is an explicit family, not a global slice theorem for arbitrary
three-dimensional metrics. At w=2 the displayed longitudinal inverse
and construction are unavailable.

Set b=beta and K(b)=<exp(-b v)>. The new log-volume coordinate is
v_w=v(phi)+log(phi')/3. Its exact spatial mean is

<v_w>_x=-(1+b/3)K'(b)/K(b)+log K(b)/3.

For any scalar N, <N_w>_x=<N exp(-b v)>/K(b).
At w=2/3, b=0 and b'=9/4. Thus

partial_w <v_w>=-(9/4)(<v^2>-<v>^2),
partial_w <N_w>=-(9/4)(<vN>-<v><N>).

These identities hold for each admissible smooth profile, independently
of a quantum state. They provide an exact nonlinear check of the
second-order coordinate-mean effect; they do not supply the interacting
source-functional Nielsen vector.
