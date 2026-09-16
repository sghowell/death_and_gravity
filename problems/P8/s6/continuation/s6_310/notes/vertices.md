# Canonical higher metric and Einstein vertices

Use eta=(+---), all momenta outgoing, lower-index metric polarizations.
Write tr(A)=tr(eta A) and g=eta+2 sum_i x_i A_i/sqrt(kappa).
Every vertex is a full mixed coefficient in distinct x_i, not a
single-field Taylor coefficient without its polarization factorial.

At unit metric coupling the determinant-density coefficients are
d0=1, d1(A)=tr(A),
d2(A,B)=tr(A)tr(B)-2tr(eta A eta B), and

d3(A,B,C)=sum_perm[
tr(A)tr(B)tr(C)/6 - tr(A)tr(eta B eta C)
+4 tr(eta A eta B eta C)/3].

The inverse metric mixed coefficient is
sum_perm (-2)^r eta A1 eta A2 ... eta Ar eta.
Convolving its subset coefficients with the determinant coefficients
gives K_r=[x1...xr]sqrt(-g)g^-1.
The scalar-pair vertex of mass-squared m2 is
-[p_lower K_r q_lower + m2 d_r]/kappa^(r/2).
The interactions -g H Phi^2/2 and C Phi^4/24 give -g d_r and C d_r,
respectively, with the same metric coupling and their matter factorials
already differentiated. The metric tensors are not TT inside graphs.

For Einstein vertices use the Gamma-Gamma density and the S304
connection/bilinear notation. At unit undoubled metric perturbation,
M1(A)=tr(A)eta/2-eta A eta and the ordered second coefficient is

M2(A,B)=eta[tr(A)tr(B)/8-tr(eta A eta B)/4]
-tr(A)eta B eta/2+eta A eta B eta.

The third connection is C3(A,B;C)=eta A eta B C1(C).
The full quartic vertex is -8/kappa times the sum over all24
permutations(A,B,C,D) of

B(M2(A,B),C1(C),C1(D))
+B(M1(A),C2(B;C),C1(D))
+B(M1(A),C1(C),C2(B;D))
+B(eta,C3(A,B;C),C1(D))
+B(eta,C1(C),C3(A,B;D))
+B(eta,C2(A;C),C2(B;D)).

The factor8 is (kappa/2)*(2/sqrt(kappa))^4; the minus sign comes from
the two Fourier derivatives. The cubic vertex remains the independently
checked S304 -4/sqrt(kappa) expression. Linear de Donder gauge contributes
no higher gauge-fixing interaction.

A separate implementation expands the literal determinant and adjugate
of eta+2tH to check every third-order inverse-density matrix entry.
Another direct indexed calculation forms Christoffel symbols and their
contraction from the metric inverse, without using the inherited
connection or bilinear helper. Inclusion-exclusion over subsets extracts
the mixed fourth coefficient. At the recorded nonzero test it gives3286,
exactly matching the new quartic vertex. These numerical tensor
calibrations supplement the displayed general coefficient derivation.
