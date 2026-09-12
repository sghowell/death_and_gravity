# Identical-state normalization and the first elastic coefficient

Use positive canonical one-particle normalization, the fully labelled invariant amplitude A, and S=1+iT. In the center-of-mass frame of two mass-one scalars,

dPhi2=beta dOmega/(32pi²), beta=sqrt(1-4/S).

The full sphere contains both permutations of two identical final particles. There is exactly one factor1/2! in the completeness sum. The forward optical identity is therefore

2 Im A(S,0) >= (1/2!) integral dPhi2 abs(A(S,x))²,
Im A(S,0) >= beta/(64pi) integral_-1^1 abs(A(S,x))² dx.

The inequality includes all other positive intermediate channels and is an equality for the leading elastic contribution. No additional final factorial is inserted. Normalize the identical incoming and outgoing two-body channels with their respective1/sqrt(2) factors. Equivalently, define even partial waves by

t_l=beta/(64pi) integral A P_l,
A=(32pi/beta) sum_even(2l+1)t_l P_l,
S_l=1+2i t_l.

This gives t0=beta Q0/(32pi), t2=beta b/(240pi). It can also be verified directly by comparing the angular norm with partial-wave unitarity. Rotational invariance separates l; exact unitarity of the full channel matrix implies abs(S_l,elastic)<=1 for its elastic diagonal entry even with inelastic states. Consequently

Im t_l >= abs(t_l)²,
(Re t_l)²+(Im t_l-1/2)² <= 1/4.

No assumption of an exactly elastic or normal scattering matrix is needed. The independent two-channel orthogonal matrix fixture verifies the strict inelastic case. These statements concern an existing exact S matrix, not the real tree by itself.

The complete first elastic absorptive coefficient determined by the original tree is

rho1=beta[Q0²+4b²/45]/(32pi)
     =(32pi/beta)(t0_tree²+5t2_tree²).

At the first four-scalar loop discontinuity the allowed two-particle tree cut is the two-scalar channel. The original vector-source vertices cannot make a two-particle four-leg production tree, and gravity is not in the nongravitational limit. Higher-field tadpoles can contribute real local terms but do not alter this elastic cut identity. This is not a calculation of the complete real one-loop amplitude or a bound on higher orders.

A constant labelled vertex c is a useful independent normalization control: Im A1=beta c²/(32pi) gives Im t0,1=(beta c/(32pi))². Applying a16pi projection to that same labelled amplitude makes the left side only half the alleged tree-wave square. The missing identical-state factor would therefore change the named matching contradiction. The subsequent dispersion argument uses the EXACT optical inequality plus a full-amplitude error premise, not a loop truncation without control.
