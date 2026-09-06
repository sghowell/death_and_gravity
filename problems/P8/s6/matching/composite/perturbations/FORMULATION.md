# P8-S6.7.COMPOSITE — a bounded rolling tensor/vector screen

This child pins and replays S6.6.COMPOSITE, certificate
`7708eebb7f10be8bad523029aa7a53921c943fd6f7df66691b71651d090ec06f`.
It does not modify that background action or its evidence. The prescribed
physical metric is still the **composite** metric, not the old CD/M1 matter
metric. The one canonical scalar is shared through that coupling; the
separate-source HR exclusions are not applied to it.

The domain is the common spatially flat, positive-square-root FLRW branch,
with smooth positive finite lapses and scale factors, positive Einstein
coefficients G,F, constant real HR beta_n, m4>0, and positive composite
couplings alpha,beta. The source rho,p is the actual homogeneous canonical
scalar source of the pinned backgrounds. Perturbations are evaluated on a
background satisfying all source-aware lapse/scale equations.

The gate establishes:

1. A literal diagonal-exponential TT quadratic action, including the
   composite matter pressure. Its relative stiffness is
   `mu=y{m4[beta1+beta2(y+c)+beta3*c*y]-alpha*beta*r*s*p}`,
   where y=b/a, c=Nf/Ng, r=alpha+beta*y and s=alpha+beta*c.
2. An independent 2x2 square-root derivation of the finite shift coefficient
   `Xi=2*y^2*Z_v/(c+y)`, with
   `Z_v=P+alpha*beta*r^2*rho+alpha*beta*r*s*y*(rho+p)/(c+y)`.
   No division by mu is used at its zero slice. P=m4(beta1+2beta2*y+beta3*y²).
3. Exact positive completed-square vector constraint elimination at any
   k>0, Xi>0. For the literal source convention q=k(E-S),
   `L=K*qdot^2-U*q^2`,
   `K=G*a^3*C/(16Ng)`, `U=Ng*a^3*mu/16`,
   `C=[1+G*c/(F*y^3)+G*k^2/(a^2*Xi)]^-1`.
   The physical high-k principal speed squared is `(r/s)^2*mu/Xi`.
4. Exact physical-clock and time-dependent canonical normalization
   identities. A zero of mu alone is not a statement that a rolling
   canonical frequency or a heavy-mode gap vanishes.
5. On the two frozen S6.6 solutions, a negative vector principal coefficient
   in some punctured physical-time neighborhood of the bounce. The vector
   inertia remains positive and regular there. The proof is analytic, with
   exact nonzero leading jets, not a numerical sampled sign test.
6. A guard against generalizing that exclusion: with independent parent
   interaction scale m and CD duration tau, the leading coefficient changes
   sign at m*tau=sqrt(24). A positive leading coefficient is not a stability
   or matching theorem, and the equality case needs higher jets.

The normalization dictionary and the direct elimination of source equation
4.8 are part of the proof. The printed equation 4.12 is not silently used
as a second normalization oracle in the absence of another specified field
transformation.

The result excludes only the formal vector-principal-health requirement for
the two frozen analytic candidates near their bounces. It does **not**
establish a large accumulated instability below a known cutoff, a numerical
width of the negative interval, a full coupled scalar action, rolling
canonical spectral gaps, nonlinear interactions, a background-dependent
cutoff, generic absence of the BD mode, positivity matching or a healthy UV
completion. It neither excludes every composite model nor completes P8.

Proof: [notes/modes.md](notes/modes.md).
Source audit: [notes/sources.md](notes/sources.md).
Report: [certificates/composite-modes.json](certificates/composite-modes.json).
