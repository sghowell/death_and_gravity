# Exact rational positive-polynomial decompositions

Each certificate targets the polynomial

    P=Fcoarse-U sum(alpha)/4.

The decomposition uses nonnegative monomials and positive
weights multiplying squares. Coefficients are exact rationals;
a numerical sample or optimizer success is not proof.

After the explicit tadpole bases below, every remaining
negative monomial has three distinct variables, -b xi xj xk.
A term w xk (xi-xj)^2 supplies a negative cross coefficient
-2w and consumes w of each available positive coefficient
of xk xi^2 and xk xj^2. A deterministic capacity allocation
uses the natural polynomial order and the three possible
choices of the outside variable. The code requires every
negative coefficient to be fully canceled, then expands the
complete identity and verifies that every residual monomial
coefficient is nonnegative. All 192 refinements succeed.

For the canonical tadpole family label the original light
parameters d (self-loop), a (external-pair edge), b,c (the
two sides incident on the self-loop vertex), and A=a+b+c.
At zero heavy parameters the necessary larger square base is

    P0=3d[(a-b-c)^2+dA]/4.

This accounts for repeated-variable negative terms that the
three-distinct-variable capacity step cannot certify.

If the self-loop vertex's pairing is 2 or 3, its light
self-loop is split into a mixed light-heavy bubble. Let h
be this particular heavy parameter. The exact base at zero
other heavy parameters becomes

    P0 + 3h(a-b-c-d)^2/4
       +63hdA/4 +63h^2(A+d)/4.

The first polynomial is recovered at h=0. Both displayed bases
are explicitly positive squares plus nonnegative monomials.
The native calculation checks the restriction identities before
allocating the remaining heavy-dependent terms.

The report records each choice tuple, tree/two-forest counts,
square allocation, positive-remainder counts, and SHA-256
digests of canonical polynomial terms for U, F and P. Full
polynomials and square bases are rebuilt from pinned source,
and the entire symbolic remainder identity is checked.
Hashes identify outputs; it is the exact algebra and positive
coefficient checks, not the hashes, that establish the bound.
