# Nonzero relational clock residue

Use S6.99's exact canonical gamma packet chart, its symplectic
normalization E/sqrt(k), and the SAME complete global equation.
Let D=diag(kappa_c,1), omega=diag(c/a,1/a), V=[1,0;-ell,1].
The action-normalized negative modes are

    S_minus=1/sqrt(2) *
       [V*(D*omega)^-1/2;
        -i*V^-T*(D*omega)^1/2].

Native identities give S_minus^dagger*i*Omega*S_minus=I and
symplectic orthogonality to the conjugate modes. They diagonalize
the actual leading canonical generator with frequencies
-i*k*c/a and -i*k/a. The generic kinetic formula is matched to
the actual global K using w=-ell*Lambda and
kappa_c=2*(J+delta_J)/Lambda^2.

The complete order-zero generator is diag(E,-E^T), with real E.
Each diagonal entry of the mode transport
S_minus^dagger*i*Omega*(L0*S_minus-S_minus') vanishes.
The native check keeps every derivative of a,c,ell,kappa_c;
it does not set those derivatives to zero. This also follows
by writing a normalized column as (q,-ip): its diagonal L0
term is p^T E q-q^T E^T p=0, and its derivative term is
(p^T q)'=0. Subleading transport is retained by the full
asymptotic construction.

Returning to the physical density coordinates, the leading
chi components have spectral weights

    W_chichi ~ hbar/(2*kappa*a^2*k) *
                   [r_c*clock_phase + matter_phase],
    r_c=ell^2/(kappa_c*c).

The inverse-square-root frequency and scale factors come from
the ACTUAL symplectic map. They are not assigned oscillator
normalizations. At u=0 one obtains r_c=1/(1215*c0).

There is a strictly positive clock weight on the whole gamma slab.
There ell>=1/[10*(17/16)^6], |Lambda|>=1/5,
J+delta_J<36+1/50 and c<1. Therefore

    r_c > {1/[100*(17/16)^12]} *
           (1/25)/[2*(36+1/50)] > 0.

All bounds are from the actual polynomial and chart domains.
The reduced matter observable has not accidentally projected
out the clock branch.

The all-orders phase construction in S6.99 is elliptic in this
nonzero scalar projection near a clock covector at coincidence.
To see the wavefront conclusion, localize to an angular cone
and to a small neighborhood of a gamma-slab point. The clock
phase and the matter phase have distinct time gradients.
The localized clock integral has a nonzero order-minus-one
leading scalar amplitude; stationary localization to its
phase covectors leaves that algebraic leading term. It cannot
decay to every order. The matter term is nonstationary in that
microlocal neighborhood, and all smoothing or lower-symbol
corrections cannot cancel the leading clock coefficient.
Thus the clock diagonal wavefront belongs to W_chichi.

Its transpose has the opposite time-frequency sign. The actual
commutator consequently also has this clock singularity. But
a clock covector has g^-1(p,p)=(1-c^2)*k^2/a^2>0, whereas
ordinary free-matter Hadamard covectors have g^-1(p,p)=0.
Any purported ordinary one-metric Hadamard covariance would
give a commutator with only that metric's null wavefront and
contradict this nonzero clock coefficient. The incompatibility
is with that replacement field/state, not with the proved
generalized Hadamard theory or with matter causality.
