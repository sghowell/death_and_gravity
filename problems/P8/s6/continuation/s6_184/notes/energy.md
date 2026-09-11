# Constrained causal response and actual mean scalar force

At physical Minkowski g, retain all three vector polarizations,
the temporal constraint and the full source contact. Write

    A=sqrt(kappa)W/m, J=sqrt(kappa)mS,
    A0=(J0-div pi)/m^2.

The unchanged S6.178 Hamiltonian specializes to a=1. Its
positive reference energy is
E0=integral[pi^2/2+Fij^2/4+m^2 A_sp^2/2+(div pi)^2/(2m^2)].
It is not the full source-driven Hamiltonian.
The exact energy identity and spatial integrations by parts give

    (sqrt(2E0))'<=F_J,
    F_J=sqrt(kappa)m F_S,
    F_S^2=||S_sp||2^2+||div S_sp||2^2/m^2+||grad S0||2^2/m^2.

The last two source terms do not vanish for a general source.
Using ||grad S||<=2C U_Phi gives
||div S_sp||<=4C U_Phi and ||grad S0||<=2C U_Phi.
Thus F_S^2<=C^2(1+20/m^2)U_Phi^2<4C^2U_Phi^2.
Regularizing sqrt(2E0) at zero energy justifies the
differential inequality with prepared zero data. It follows that

    E0(t)<=2kappa m^2 C^2 (integral_past U_Phi)^2.

The actual FREE coherent vector energy instead contains
(J0-div pi)^2/(2m^2). Hence
Efree<=2E0+||J0||2^2/m^2, giving

    Efree(t)<=4kappa m^2 C^2 (integral_past U_Phi)^2
              +kappa C^2 U_Phi(t)^2.

This is still not the full source-interaction stress.
At the anchor the three coefficients above are respectively
below1e-2376, 2e-2376 and4e-2383.

Energy controls W_sp, and the temporal constraint controls W0.
The Euclidean four-component norm obeys
||Wbar||2<=2m integral_past F_S+||S0||2
          <=4m C integral_past U_Phi+C U_Phi(t).
Thus ||S-Wbar||2<=4m C integral_past U_Phi+2C U_Phi(t).

The expectation of the UNINTEGRATED scalar source-sector
variation is exactly
Force_Phi[eta]=kappa integral (S-Wbar).DS_eta.
The fixed-metric determinant is Phi-independent and the
centered linear Gaussian mean vanishes. Evaluate Wbar at
the causal solution; do not vary a single-branch S Gret S.

Cauchy-Schwarz and interval length<=1 yield
|Force_Phi[eta]|<=kappa C^2(4m+2)||Phi||J3||eta||J3.
At the anchor this coefficient is below1e-2378.
Every energy/force coefficient scales as kappa0/kappa
along the actual fixed-canonical family.

No real-time pole denominator is inverted in an L2 norm.
The finite-time energy argument keeps every spatial momentum
and all temporal resonances. It is not a whole-time response
bound, on-shell scattering remainder or physical cutoff.
