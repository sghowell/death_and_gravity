# Literal complete current and exact polynomial arithmetic

S315's canonical metric action fixes all vertices; only EH3 and EH4
are needed for three leaves. The new polynomial implementation expands
the same inverse metric, density and Christoffel contractions. For each
complete source it sums ALL labeled set partitions. The root count is4,
and every pair count is1. Every source Ward identity is checked as a
polynomial identity before the inverse is applied.

Ten symmetric root variables collect the complete source in one vertex
evaluation. Each surviving monomial must be linear in exactly one root
variable. Off-diagonal coefficients are divided by2. The independent
literal action evaluator checks EH3/EH4 and complete pair/root sources.

For source J and total Q, put
 S=-eta J eta+eta tr(eta J)/2.
The spatial temporal inverse is exactly
 [Q0^2 Sij+Qi Q0 S0j+Qj Q0 Si0+Qi Qj S00]/(Q^2 Q0^2).
Its time row and column vanish. A generic symmetric source calculation
checks all16 components against the original inverse then projection;
this operator identity does not assume conservation.

Stereographic leaves have rational momenta and TT tensors. Multiplying
all momenta by an angle-only common denominator leaves the current
unchanged: each two-derivative vertex is paired with one inverse.
The implementation retains exact known denominator factors. It may
leave removable factors in place; no reduced-form assumption is made.

Known-factor division uses unchanged SymPy polynomial division, moving
one variable to the outer ring only when its leading coefficient is a
nonzero rational unit. Every quotient satisfies P=q*f+remainder in the
original ring. Factored denominator assembly also checks q*Dbranch=D.
Multiplication by ab/(a+b+1) checks an independent denominator identity.
These are representation choices, not global backend replacement.

The production code reconstructs every polynomial; it imports no
private cache. All private/direct/native/ordinary/CLI paths use original
SymPy. The exact-GCD adapter is restricted to the complete FULL runner.
