# Actual canonical Proca action, constraint and physical readouts

The retained action is unchanged:

    kappa sqrt(-g)[-zeta F(W)^2/4+(W-S)^2/2].

Fix zeta=10^-6, allowed by S6.174, and kappa=10^800. Define

    A=sqrt(kappa*zeta)W,  m^2=1/zeta=10^6,
    J=sqrt(kappa/zeta)S.

The canonical action is ordinary Proca plus -A.J+kappa S.S/2.
The source normalization and contact term are checked exactly.
The chosen canonical three-mode measure is stated in renormalization.md;
this does not silently fix all other parent-field Jacobians.

For a general positive homogeneous lapse N and scale a, at fixed comoving
momentum k, the coordinate-vector transverse action is

    L_T=a A_dot^2/(2N)-N a(m^2+q)A^2/2, q=k^2/a^2.

For the longitudinal spatial amplitude, retain its temporal component.
One real sine/cosine Fourier-pair convention gives

    L_L=a(A_dot-k A0)^2/(2N)+a^3 m^2 A0^2/(2N)
           -N a m^2 A^2/2.
    A0=k A_dot/(k^2+a^2 m^2).

Its exact reduced kinetic coefficient is g_L^2=a m^2/(q+m^2);
g_T^2=a. Both have reduced Hamiltonian

    H=N[pi_A^2/(2g^2)+g^2(q+m^2)A^2/2].

Vary N and a while holding coordinate A,pi_A and comoving k fixed,
BEFORE making the metric-dependent canonical substitution
v=gA,p=pi_A/g. Its fixed-metric symplectic determinant is one.
The physical mode energy and isotropic pressure are

    rho_T=rho_L=[p^2+(q+m^2)v^2]/(2a^3),
    P_T=[p^2+(q-m^2)v^2]/(6a^3),
    P_L=[(1+2z)p^2-(q+m^2)v^2]/(6a^3),
    z=q/(q+m^2).

These formulas include two transverse modes and one longitudinal mode
when summed. A pressure coefficient need not be positive individually.
The expressions are derived from literal metric variations and compared
to the original ordinary readout matrices, not inferred from equal
clock covariance. At k=0 all three canonical mode equations agree;
the isotropic pressure is understood after the polarization sum.

At N=1 the canonical evolution is v'=p+d v,p'=-omega^2 v-d p,
with omega^2=q+m^2, d_T=H/2 and d_L=H(1/2+z). Hence

    v''+[omega^2-d'-d^2]v=0.

For a=(1+t^2)^2 these are exactly the frozen S6.55 clock equations.
The calculation retains z'=-2Hz(1-z) at fixed k. It does not replace
the longitudinal equation by a transverse oscillator.

The old -+++ metric is minus the P8 +--- metric. With identical lower
A components, their Levi-Civita connections and volume elements agree;
the Maxwell double contraction is unchanged and the mass-term sign
changes with the metric. Thus the two physical ordinary actions,
coordinate-mode Hamiltonians and positive canonical covariances agree.
This is not a four-scalar positive-fiber quantization.

The same positive-frequency mode convention is exp(-i*omega*(t-s)).
The sign of a raised wavefront covector changes under the overall metric
signature reversal; the old theorem's future-covector label is mapped
with that reversal, not copied while changing its definition.
