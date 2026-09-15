# Exact corrected finite-model formulation

## Original boundary and units

The original S220 scalar Lagrangian is per kappa a_bar^3 in the
reference density variables v,sigma,pv,ps. Its integration by parts
removes the generating function

    f = -9 Hbar v^2 + 3 ellbar v sigma.

Thus the RAW original ADM density momenta are
r=pv-18Hbar v+3ellbar sigma and m=ps+3ellbar v. The actual reference is

    a_bar=(1+u^2)^2,
    Hbar=4u/(1+u^2),
    ellbar=1/(10 a_bar^3).

With unit-CCR configuration Q=sqrt(kappa)(v,sigma) and prepared momentum
P=sqrt(kappa)a_bar^3(pv,ps), the raw canonical map is

    z_raw = Sb(u) z_prepared,
    Sb = [I,0; C,I],
    C = a_bar^3[-18Hbar,3ellbar;3ellbar,0],
    F_b = Q^T C Q/2.

The off-diagonal C entries are exactly3/10. The physical finite field
generator is kappa a_bar^3 integral f. Both descriptions agree because
the six real Fourier functions are orthonormal. The SAME block acts on
each scalar pair; tensor, heavy and all Proca channels are unchanged.

This is an explicitly stipulated exact LINEAR canonical extension on
the full finite reduced chart. The raw source is nonlinear and is
pulled back exactly. No unproved full nonlinear boundary generating
function or nonlinear metaplectic covariance is claimed.

## Corrected state, time generator and full Hamiltonian

Use the original prepared Gaussian psi0, its full initial whitening S0,
prepared reference flow Mref and exact propagator Uref. Set

    Ub(u)=exp(i F_b),
    psi_raw(0)=Ub(0) psi0,
    Uref_raw(u,0)=Ub(u) Uref(u,0) Ub(0)^-1.

The raw reference flow is Sb(u)Mref(u,0)Sb(0)^-1 and its initial whitening
is Sb(0)S0. This transports the original stipulated state; no vacuum is
selected again, cross covariance discarded or Gaussian substituted.

Let hraw(u,Y,zraw) be the complete original absolute finite reduced
Hamiltonian kappa integral[V Hbar(Nstar,Z)], with every fixed source,
primitive, temporal/spatial root, shape/cotangent inverse, Gauss
contact and generated harmonic retained. In prepared coordinates use

    hcor = hraw(u,Y,Sb(u)z) + partial_u F_b(z).

The sign follows both from the full canonical one-form and the exact
Schrodinger conjugation. The added quadratic energy is

    partial_u F_b = -9 kappa a_bar^3(Hbar_dot+3Hbar^2) integral v^2.

The cross-coefficient derivative vanishes, not the v^2 term: at u0
its unit-CCR value is -36 Q_v^2. The generating function is fixed
reference-time data, independent of LIVE homogeneous Y.

For fixed original coreR=1e20 and each original c1,c2 cutoff define

    g=(hcor-href)(u,Y,Mref z),
    g0=hraw(u,Y,0),
    b_c=chi_c(g-g0),
    K_C=Q_V[(1-D)b_c],  K_W=OpW(b_c),  D=Delta_w/4.

This changes the old raw/prepared physical model, not the definitions
of the two orderings. All derivatives of the COMPLETE localized b_c,
including the boundary term, are retained. Neither polynomial heat
truncation nor a symbol-supremum Weyl-norm replacement is used.

## Coupled homogeneous equations and volume

The live classical density coordinates are
Y=(alpha,p,pm,eta=1e100 h,ph); M1 is reconstructed by its full cyclic
rate. The absolute homogeneous canonical one-form is

    kappa Vol[PV d alpha + PM d M1 + PH d eta/1e100],
    PV=3a^3p, PM=a^3pm, PH=a^3ph.

Keep the entire original scalar-center classical force plus every
canonical quantum-force derivative of the corrected K. No live force
is removed with g0. The configuration-dependent Ub is not a scalar
phase. The physical state is

    psi_raw(u)=Ub(u)Uref(u)exp[-i integral g0(Y)]phi(u),
    i phi'=K_Y phi.

All three homogeneous spatial-vector and five traceless-shape
coordinate/momentum pairs satisfy their equations at zero by actual
cubic symmetry, not a two-TT-zero-mode truncation. The heavy scalar
remains live. The fixed profiles are not replaced by evolving means.

Pull the complete volume symbol through Sb and Mref:

    F_Y=average exp(3v)R_full(u,Nstar)^(-3/4),
    F0(Y)=R_full(u,N0(Y))^(-3/4),
    F_ext=F0+chi_c(F_Y-F0).

In each ordering, nu=exp(3alpha) times the expectation of F_ext.
State, physical readout and outside-core positive coherent effect
are conjugated by the SAME Ub and Uref. No generating-energy term is
added to the volume.

## Evaluated domain and exact claim

Keep L1,P1e64,kappa1e800,zeta1e-6 and all48 nonzero-mode pairs.
T=1e-180, real homogeneous radius1e-390, complex radius1e-122,
Cauchy radius1e-124, source box1e-120 and lapse ball1e-115.
Rstar1e150 is only a larger complex ANALYSIS radius.

The two corrected raw momentum rows are doubled; every other row is
unchanged. Complete spatial/source bounds and the boundary integral
give centered amplitudes AH1e172,AF1e-510. All206 phase derivatives and
first/second homogeneous derivatives enter exact operator estimates.

Each corrected finite hybrid has a unique coupled path on[-T,T].
The volume is C1, has relative error<1e-380, both endpoint gaps>5T^2
and all minima inside |u|<1e-188. Positive core leakage is<1e-6.
The complete two-row coupled comparison gives the numerical bounds
in README.md. State-vector comparisons factor each OWN global g0 phase
only; they do not compare differently phased un-factored states.

Original V/G/B/P8 remain OPEN. This finite physical bridge is repaired
without claiming a continuum theory, homogeneous quantum state, unique
strict minimum, physical UV matching or global completion.
