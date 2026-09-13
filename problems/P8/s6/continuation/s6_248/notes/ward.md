# Both Gaussian Ward legs, the full clock contact and the force coordinates

All nonlocal scalar response is moved to synchronous spatial directions by the complete metric Ward identities, not by discarding lapse/shift equations. For opposite Fourier signs, with xi=(eta,0,0,-sign i c/P),

D_phys=(eta', (0,0,-sign i(c'+q eta)/P), 2(w+H eta)I),
D_syn=(0,0,2wI-2cPi), Pi=diag(0,0,1).

The literal first ADM metric jets verify D_phys=D_syn+D_gauge on both detector and source legs. The separate 1/P in xi is not bounded before combining the full expressions.

For each actual species, E is its nonzero weight-one contravariant Gaussian current,
E=diag(-a³rho/2,-aP/2,-aP/2,-aP/2).
With h the first metric jet and q2 the full mixed second metric jet, the ordered local completion is

W_s(D,xi_G;E)=h(D):Lie_density(xi_G)E+E:q2(D,G_gauge),
W_d(xi_D,G_syn;E)=-E:Lie_covariant(xi_D)h(G_syn)
                              +E:q2(D_gauge,G_syn).

Apply source Ward with the FULL detector first, then detector Ward with the SYNCHRONOUS source. Add the distinct Gaussian common-clock contact
3a³P(4delta²-3delta)eta_D'eta_G'
exactly once. The S246 scalar kinematics and canonical density identity are used directly, with the full actual Proca-plus-H mean substituted only after the generic tensor-density calculation. Linearity proves that the two species' local completions sum exactly. No vector norm constant is transferred to the scalar sector.

The complete expression has fourteen bilinear jet monomials, at most one derivative on each time leg, degree at most two in P, no uncancelled Fourier phase and no denominator containing P. Thus its Euler contribution is local of order at most two, uniformly on each fixed external ball including the small nonzero-P limit. The first response does not need spatial derivatives of a selected state or a regular polarization chart at a zero internal leg.

The fixed profiles do not belong in the metric-only Gaussian Ward identity. Remove the full spatial heavy profile before the heavy Gaussian reconstruction and restore its full clock profile only in the coefficient action, as in assembly.md. Gaussian and profile nonlinear clock contacts are opposites at the full reference, but both are kept until the whole equation is assembled. Cancellation of the total reference mean never makes the Gaussian one-current individually zero.

For force recovery use the physical scalar amplitudes s=(n,zeta,b), whereas the coefficient action uses x_clock=(n,v,b). Their first-jet map is
A(t)=[[1,0,0],[delta(t),1,0],[0,0,1]], s=A x_clock.
The scalar pure-heavy clock response is
R_H,clock(t,s)=A(t)^T Rhat_H(t,s) A(s),
where Rhat_H includes the distinct clock nn contact. Therefore

Rhat_H(t,s)=A(t)^(-T)[R_S246,clock-P_H,clock](t,s)A(s)^(-1),
Qbar(t,s)=[kappa a(t)³]^(-1)[Rhat_Proca+Rhat_H](t,s).

Every detector/source matrix and the OUTPUT density retain this order. A same-time clock matrix commutes with scalar same-time density multiplication; it does not justify moving the density through a two-time response. Independent unequal-time matrix fixtures detect a missing leg and a density evaluated at the source time.

Start the tensor identities on smooth Fourier data away from P0. The combined expressions and finite-ball bounds extend in scalar-amplitude Sobolev norms by density. Multiplication by1/P maps L2 on a bounded three-dimensional ball locally into L1, so the physical shift is a tempered distribution; no identical infrared H^r bound for the shift vector is asserted. The literal homogeneous P0 constraint system is distinct from this almost-everywhere Fourier-amplitude extension.
