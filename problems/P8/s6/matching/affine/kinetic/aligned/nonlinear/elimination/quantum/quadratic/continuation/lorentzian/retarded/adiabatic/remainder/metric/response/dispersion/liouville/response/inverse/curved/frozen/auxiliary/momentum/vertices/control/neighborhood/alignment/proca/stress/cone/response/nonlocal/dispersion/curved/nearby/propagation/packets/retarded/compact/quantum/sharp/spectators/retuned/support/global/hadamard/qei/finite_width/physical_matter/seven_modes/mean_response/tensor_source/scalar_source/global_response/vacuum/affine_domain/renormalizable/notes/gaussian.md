# Exact finite-regulator integration retains mixed loops

At a fixed finite Euclidean regulator use positive symmetric
K>=mu^2 I, K_m>=I, the light vector phi and J_i=phi_i^2.
The regulated action is

    S=phi.K_m.phi/2+H.K.H/2+G H.J/2
      +lambda4 sum_i phi_i^4/24.

With H_*=-G K^-1 J/2, completing the Gaussian gives

    S=(H-H_*).K.(H-H_*)/2+S_eff,
    S_eff=phi.K_m.phi/2+lambda4 sum_i phi_i^4/24
          -G^2 J.K^-1.J/8.

The exact Gaussian determinant is independent of phi.
Because K^-1<=mu^-2 I, S_eff>=phi.K_m.phi/2+
delta sum_i phi_i^4/24. Both finite-dimensional partition
functions converge, and normalized light correlations agree.
A finite-volume or continuum limit is not asserted by this
finite-regulator fact. Bare and renormalized positivity must
also not be silently identified across counterterms.

Differentiate the full action BEFORE imposing H=H_*.
Its Hessian blocks are A=S_phiphi, B=G diag(phi), K.
Differentiating the stationary action gives exactly

    S_eff,phiphi = A(H_*)-B K^-1 B^T,
    det S_full''(H_*) = det K det S_eff''.

The native two-site calculation uses independent generic
entries of both kinetic matrices, all fields and both
couplings. It verifies the original stationarity equations,
stationary action, full Hessian blocks and Schur complement.
An independent 24-permutation four-by-four determinant
checks the literal determinant identity without a symbolic
matrix-pivot prescription.

The nonlocal inverse K^-1 is retained inside every light
loop. Mixed heavy/light loops are present in the resulting
light functional even though det K is field-independent.
Replacing K^-1 by a finite low-momentum polynomial before
an unbounded loop integral is not this exact identity.
The one-loop full-Hessian reduction is also the standard
resolved determinant construction reviewed in
[literature.md](literature.md).

A constant heavy tadpole, a constant heavy wavefunction
normalization and renormalized heavy mass keep the same
Gaussian form, with a shifted source and kinetic kernel.
The power-counting proof lists which divergent structures
are actually required. No field-dependent quantum-gravity
measure or old affine measure is identified with this one.
