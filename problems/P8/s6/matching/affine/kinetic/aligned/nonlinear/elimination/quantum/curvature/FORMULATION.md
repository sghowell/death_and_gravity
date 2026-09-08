# S6.49: retained-vector local curvature counterterms

Certificate-gated limited component; original P8 OPEN. Do not alter S6.42 or any
frozen ancestor. On the exact clock the fixed-light vector Hessian
is ordinary constant-mass Proca on the physical rolling metric.
Derive its local covariant pole, using a stated positive Euclidean
Hodge-Laplacian prescription and epsilon_DR=(4-d)/2, then evaluate
the invariant coefficient functions on the actual Lorentzian FLRW
geometry. Keep total derivatives and their boundary qualification.

The on-clock Hessian does not determine off-clock scalar/metric
variations: a_N,b_N and higher mass derivatives do not vanish.
Accordingly a background counterterm evaluation is not the full
quantum stress, scalar constraint Jacobian or renormalized bounce.
No finite curved or in-in remainder is asserted from heat-kernel
asymptotics alone. The original V/G/B acceptance gates remain open.

Required controls: the full vector-minus-scalar trace, exact tensor
contractions, mass-expansion signs, flat S6.47 pole matching, and
an independent sphere spectral expansion with excluded modes explicit.
Use the actual metric, not a replacement de Sitter solution. Any
bound is on the displayed local coefficients, not the full loop.

## Exact result and domains

Use D1=-nabla^2+Ricci on one-forms and D0=-nabla^2 on scalars.
The covariant dimensional prescription gives
Gamma_1loop=1/2 Tr log(D1+m0^2)-1/2 Tr log(D0+m0^2), with
scaleless contact traces set to zero only in this prescription.
The local heat coefficients of their difference, excluding (4*pi*s)^-2,
are a0=3, a2=-R/2 and

    a4=-R^2/8+29*Ricci^2/60-Riemann^2/15-box(R)/15.

The Euclidean pole is -[3*m0^4/2+m0^2*R/2+a4]/(32*pi^2*epsilon_DR).
Analytically continued local invariant counterterm functions are
evaluated on the original physical Lorentzian metric, not interpreted
as a finite Lorentzian determinant. With u=t/time_scale,

    time_scale^4*a4 = -8*(77*u^4+14*u^2-3)/(1+u^2)^4,
    0 < time_scale^2*R <= 49,
    |time_scale^4*a4| <= 308/3

for all real u. The fourth-order coefficient includes its displayed
total divergence. Dropping that divergence changes the coefficient;
no vanishing boundary flux at infinite cosmological time is assumed.
For Rm=m0*time_scale>0 the absolute curvature correction to the pole
weight, divided by its flat weight 3*m0^4/2, is bounded by
49/(3*Rm^2)+616/(9*Rm^4). This is a pole-coefficient ratio, not
a bound on a finite quantum correction.

The independent unit-S4 spectral trace Tr(exp(-s*D1))-Tr(exp(-s*D0))
equals the coexact one-form trace minus the scalar constant mode.
For 0<s<=1, the difference from 1/(2*s^2)-1/s-11/30 is bounded
in absolute value by (1609/1728)*s using the sixth-order
Euler-Maclaurin remainder. This verifies the local coefficients on
that diagnostic geometry; it does not transfer the finite bound to FLRW.

Pin and fully rebuild S6.48 and its ancestry. Exact tensor identities,
continuous coefficient and spectral-error proofs, strict source hashes,
read-only replay and rejected inputs/mutations gate this limited result.
No proof-assistant formalization or V/G/B completion is asserted.
