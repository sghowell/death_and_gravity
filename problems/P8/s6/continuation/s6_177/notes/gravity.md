# Exact Einstein boundary, normalized tensors and the full classical limit

Choose the canonical tensor convention g=eta+2h/sqrt(kappa). This is a
fixed factor2 relative to using sqrt(kappa)(g-eta); it is not a new
kappa-dependent interaction. The physical curvature and EH signs are
the original P8 +--- convention, with action -kappa R/2.

Keep the exact Einstein divergence before the limit. Its linear term is

    -sqrt(kappa) partial_mu[
        partial_nu h^(mu nu)-partial^mu trace(h)].

The full nonlinear divergence is retained as a boundary as well.
It integrates to zero for the compact/Schwartz variations used here;
no infinite-time boundary flux of the bounce is silently erased.
After this exact boundary, the EH density is

    -kappa sqrt(-g) g^(mu nu)/2
       [Gamma^rho_(mu sigma) Gamma^sigma_(nu rho)
        -Gamma^rho_(mu nu) Gamma^sigma_(rho sigma)].

Writing epsilon=kappa^(-1/2), Gamma=2epsilon Gamma1+O(epsilon^2)
gives the finite quadratic density computed directly by gravity.py.
Its Fourier-symbol form is the Fierz-Pauli expression

    [k^2 h_mn h^mn-2(k.h)_m(k.h)^m
       +2(k.h).k trace(h)-k^2 trace(h)^2]/2.

A de Donder gauge-fixing term gives the propagator iP/(p^2+i0),
P=(eta_mr eta_ns+eta_ms eta_nr-eta_mn eta_rs)/2.
Pure linear gauge directions have zero quadratic symbol.
For independent real TT time and z derivatives the density is

    h_plus_dot^2+h_cross_dot^2-h_plus_z^2-h_cross_z^2.

Thus sqrt(2)h_plus and sqrt(2)h_cross have the standard positive
half-normalized kinetic terms. The exact cosmological boundary check
independently gives the reduced EH density -3kappa a adot^2/N,
including the lapse derivative in the unreduced curvature.
Dense FULL inverse-metric/Christoffel evaluations check the approach
to this quadratic limit, not just a substituted quadratic formula.

For the remaining full action, take any compact set of canonical
C^(j+2) jets strictly inside the fixed scalar coefficient strip.
At sufficiently large kappa the metric and its inverse lie in a
compact nonsingular Lorentzian neighborhood of eta, and Y stays
inside that strip with a positive margin. The exact inverse identity

    g^(-1)-eta^(-1)=-2epsilon eta^(-1)h g^(-1)

and its differentiated versions bound all finite jet differences by
constants times epsilon. The determinant square root is analytic
on this compact metric neighborhood.

Each EH connection is O(epsilon), so after the exact boundary its
quadratic term is finite and every additional h factor is suppressed
by at least epsilon. The nonminimal -r R/2 has no prefactor kappa;
R=O(epsilon) on these fixed jets, so it vanishes. Alternatively,
integration by parts keeps its derivative-r terms explicitly and
gives the same rate. Covariant scalar Hessians approach flat Hessians,
and all metric changes in the complete fixed f and a3 terms vanish
at this rate. The two dependent Ia terms vanish as1/kappa.
The complete vector-source estimates are in source.md.

With Psi=sqrt(kappa)chi, original M1 becomes a fixed canonical
massless spectator. With A=sqrt(kappa*zeta)W, the source-free vector
action is fixed mass1000 Proca. Therefore the full classical limit,
modulo the retained Einstein divergence, is

    free Fierz-Pauli + f(Phi,Y)+a3(Phi,Y)(L3-L4)
      + (partial Psi)^2/2 -F(A)^2/4 +A^2/(2zeta).

Taylor's integral remainder on the compact jet set supplies a finite
C_j, independent of kappa, for the O(kappa^(-1/2)) remainder and its
first j jet derivatives. These are full-function bounds, not a
Taylor truncation in Phi. They need NOT be small at kappa0: the
canonical bounce jets are enormous, and the constants can depend
strongly on the chosen compact set and on kappa0. No finite-base
background, solution, physical cutoff or quantum-continuum estimate
is obtained by setting this asymptotic remainder to zero.
