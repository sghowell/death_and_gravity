# Two angular sectors, reflection, and arbitrary complex TT leaves

Rotate one of three distinct directions to the north pole and rotate
another into the positive x-z half-plane. The second stereographic
coordinate is (r,0), r>0; take the third as (t,u), with u>0 after a
possible spatial reflection. For a generic triple one can choose an
anchor with no antipodal partner. If a pair is antipodal, the third
direction is such an anchor. Great-circle and chart-boundary cases
follow by continuity wherever the original internal denominators
remain nonzero.

Two independent charts are reconstructed:
 sector0: (0,0),(r,0),(-v,u);
 sector1: (0,0),(r+v,0),(r,u);
with r,v,u>0. They represent t<0 and0<t<r, respectively.

For the third sector t=r+v>r, set
 sin=2r/(1+r^2), cos=(1-r^2)/(1+r^2),
 O=[[-cos,0,sin],[0,1,0],[sin,0,cos]],
 D=(1+rt)^2+r^2*u^2,
 v'=((t-r)*(1+rt)+r*u^2)/D,
 u'=u*(1+r^2)/D.
All of D,v',u' are positive. The exact map sends n(0,0) to n(r,0),
n(r,0) to n(0,0), and n(t,u) to n(-v',u'). Swapping the first two
leaves therefore reduces the third sector to sector0.

This O is orthogonal and energy-independent; diag(1,O) preserves eta
and the fixed temporal reference vector. The source uses only metric
contractions and momenta, with no parity-odd epsilon tensor. Complete
partitions are invariant under leaf relabeling. Induction through
the literal vertices, inverses and temporal projection proves covariance
at every required order. The symmetric pair-energy scale S is preserved.
Sixteen independent full source/current/Ward oracles check the map.

The stereographic frame U,V is real orthonormal. Plus=UU-VV and
cross=UV+VU are Frobenius orthogonal with squared norm2.
For A=alpha*plus+beta*cross and ||A||=1,
 2(|alpha|^2+|beta|^2)=1, hence |alpha|+|beta|<=1.
Trilinearity extends the maximum of the eight computed basis bounds
to arbitrary complex unit TT leaves. The bound is on Frobenius norm,
which is O(3)-invariant; an entrywise l1 majorant is used only inside
each chart. No common global polarization frame is assumed.

Exactly coincident-pair propagator poles remain excluded. Uniform
approach bounds do not assign a physical value at those points.
