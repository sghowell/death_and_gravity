# Precisely one unknown constant is replaced

S6.180 derived the actual current-parent retained prepared
homogeneous response, without the old stress-canceling profile:

    I4(T+gamma Q)=diag(A(t),gamma F_m(partial_t^2))+V,
    B0=diag(A^-1 multiplication,gamma^-1 K_m convolution),
    |A^-1|<=a0=15625/6144, gamma=1/(64pi^2 kappa).

Its full curved remainder has a finite but unevaluated row
majorant C(1+|log r|). Previously both C and ||K_m||L1 were
unevaluated. The new result permits the explicit substitution

    Kbar=300000000000,
    beta=C(a0+Kbar/gamma),
    lambda=(4beta+1)^2.

The same weighted Volterra argument gives contraction below
2beta/(4beta+1)<1/2 and adapted C0 bound

    ||z||<=2 exp(lambda) max(a0,Kbar/gamma)||I4 g_adapted||.

This choice is valid because the actual kernel norm is below
Kbar. The equality of the substituted formulas is checked
against S6.180. This is not a new matching prescription.

C has NOT been numerically bounded. Nor is exp(lambda) shown
small, useful for nonlinear backreaction, or an optimal inverse
bound. In fact this construction is potentially extremely
large. The previously established local C1<15000 for physical
lapse recovery is unchanged and remains a different constant.

No background residual is passed to this prepared inverse:
the reference quantum one-point stress is nonzero near the
initial neighborhood. A compatible quantum preparation and
controlled nonlinear background still need construction.
