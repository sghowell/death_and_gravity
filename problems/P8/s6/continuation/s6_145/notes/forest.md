# Both overlapping Yukawa subtractions, not a reused quartic bound

For the two self-energy words, insert the complete S6.140 scalar/gauge
proper MS fermion kernel into each marked propagator of the first
Phi two-point bubble. The mark replaces S by -S Sigma_MS S.
The S6.140 finite dimensional mass and kinetic terms remain in
Sigma_MS; it is not replaced by a local anchor.

For the vertex word write D2(q,p) for the ordered product of two
fermion propagators and K0(l)=B_b(l)S(l)^2, with b=1 for scalar
exchange and b=0 for gauge exchange. Gauge caps and their sign
are retained as in the S6.142 proper vertex. The full raw trace
G has two opposite ultraviolet-divergent Yukawa subgraphs.

Subtract BOTH K0(l)D2(q,p) and K0(q)D2(l,p). Cyclic ordering/caps
are preserved; the matrix notation in the code is a noncommuting
identity. Add back each complete regulated zero-momentum proper
vertex paired with its own MS counterterm. Their finite values
are the two local S6.139 anchors

    ups=[Y(J0+2R)-6a C_F]/Q,
    |ups|<=(2Y+6a C_F)/Q<=2(Y+4a C_F)/Q.

This is the proper vertex anchor, not a canonical Yukawa conversion.
Both finite anchors occur once. There is no product of overlapping
counterterms. The separate whole-fermion-cycle and overall local
subtractions discussed in the catalog vanish after the projection.

For y much larger than m^2+x, integrate only the paired
G-K0(l)D2(q,p), together with the separately bounded complementary
subtraction -K0(q)D2(l,p). The opposite half-domain is symmetric.
Dropping the complementary subtraction would miss a proper forest.
Integrating the raw high region independently would be divergent.

All these identities are made in the regulator first. Finite
epsilon times pole terms in the proper kernels/anchors are retained.
Remaining local overall polynomials, including their finite
dimensional terms, are killed algebraically by the on-shell map.
This argument does not calculate their coefficients.
