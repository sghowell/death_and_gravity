# Full dimensional finite local tensor

Keep d=3-2epsilon, fixed four-component trace, full angular measure,
mu=m and Q=16pi^2. At each fixed real time the exact mode matrices
have the same dimensional continuation after the a^(d/2) Cauchy
normalization. Spatial isotropy gives1/d, not an early1/3.

Let K=Mdot+HM and ell=log(M^2/mu^2). Define

    A=K^2/4-3HMK/2,
    B=-M[Mddot+(Hdot-H^2)M]/2,
    C=3M^2 K Mdot/2.
    J0=4/Q[1/epsilon-ell-2/3],
    J1=4/Q[1/epsilon-ell-16/15],
    J2=8/(5QM^2).

J0, J1 and J2 integrate respectively q^2/omega^5, q^4/omega^7
and q^2/omega^7 with the full dimensional measure. The displayed
finite expansions follow from the exact Gamma integrals, with the
MS-bar factor included. Thus

    rho_bare2=K^2[1/epsilon-ell-2/3]/Q,
    P_bare2=(A J1+B J0+C J2)/d.

The fixed-comoving derivative qdot=-Hq is essential to this pressure
coefficient. The pole action is

    [M^4-(partial M)^2-M^2 R/6]/(Q epsilon).

The zeroth-order terms and vacuum/saturated-mass references are the
same covariant potential as S6.168: its pressure is minus its energy.
For the two-derivative terms, full d-dimensional metric variation gives

    rho_CT2=[-Mdot^2-d(d-1)M^2 H^2/6-2d HM Mdot/3]/(Q epsilon),
    P_CT2=[-Mdot^2/3+2M Mddot/3
           +M^2((d-1)Hdot/3+d(d-1)H^2/6)
           +2(d-1)HM Mdot/3]/(Q epsilon).

Both poles cancel. The epsilon dependence of the counterterms adds
finite terms that a four-dimensional-first variation would lose.
The finite results are the metric variations of one common action

    L2=-(partial M)^2(ell+2/3)/Q-M^2(ell-1)R/(6Q).

In particular

    rho2_MS=[-K^2 ell-2Mdot^2/3+M^2 H^2]/Q,
    P2_MS=ell[-Mdot^2/3+4HM Mdot/3+H^2 M^2
                    +2M Mddot/3+2Hdot M^2/3]/Q
             +[2Mdot^2/3-M^2 H^2-2M^2 Hdot/3]/Q.

For L=F(t)R in the P8 curvature convention, these variations are

    rho=6FH^2+6H Fdot,
    P=-2F(2Hdot+3H^2)-2Fddot-4H Fdot.

The independent lapse/log-scale-factor calculation retains general d
before taking its limit. The flat limit exactly recovers S6.168 energy
and S6.169 pressure, including the nonzero flat curvature improvement.

After rho0+rho2 and P0+P2 subtraction the full mode remainder is
integrable near epsilon=0 at each fixed time. The compact state bound
and p^-5 complete tail give a radial O(p^-2) dominating remainder;
the projector remainders decay faster. A neighborhood |epsilon|<1/4
suffices. No uniform-regulator infinite-time exchange is asserted.
Restore local counterterms, including the explicit geometric choice
in curvature.md, before using the uniform four-dimensional bounds.
