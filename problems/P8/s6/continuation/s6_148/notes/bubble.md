# Keep both physical masses in the outer bubble

The older paired row contains exactly one heavy line with mass
squared M and one inserted light line, whose fixed pole has mass
squared one. At s=0 let

    F_D(M)=Integral_T^infinity rho_D(v) B_D(v,M;0)/(v-1)^2 dv.

Using the same normalization H and z as tadpole.md,

    F_D(M)/(C/Q)=Gamma(e)H(e) Integral_0^1
       z^(2e-1)(1-z)^(3/2-e) A(z/T) J(Mz/T,e) dz,
    J(k,e)=[1-k^(1-e)]/[(1-k)(1-e)].

The two ratios z/T and Mz/T must not be identified. The first
comes from the physical inner on-shell subtraction; the second
comes from the actual heavy propagator. At M=1 this reduces to
S6.143, but the present M is much larger than one.

Separate the same F0(e) by replacing A with one and J with
1/(1-e). The difference compact integral is holomorphic near
e=0, but multiplying by Gamma(e)H leaves a nonzero simple pole.
At e=0 its weight is w h(z/T), independent of M, where

    w=(1-z)^(3/2)/z, h=A-1, L(k)=-klog(k)/(1-k).

The complete finite difference is

    delta_F=Integral_0^1 w {
       [3-4log2+2logz-log(1-z)]h(z/T)
       -A(z/T)L(Mz/T)} dz.

The finite constant 2-4log2 in Gamma(e)H and the first
e derivative of the compact integral are both retained.
Consequently F_MS(M)/(C/Q)=47/18+pi^2/12+delta_F.
The MS pole is subtracted; its finite product is not discarded.
For 1<=M<=T/16, remainders.md bounds this whole correction.

For real 0<=s<=1, Delta=xv+(1-x)M-x(1-x)s
is bounded below by x(v-4)+(1-x)M. The gap is
x(1-x)(4-s)+4x^2. With w(v)<=4C/v and v-4>=v/2,
the convergent derivative obeys |F_MS'(M;s)|<=4C/(QT).
No local mass subtraction changes this derivative. This is the
same conservative constant used for the earlier on-shell slope,
now stated on the entire interval needed to reach s=1.
