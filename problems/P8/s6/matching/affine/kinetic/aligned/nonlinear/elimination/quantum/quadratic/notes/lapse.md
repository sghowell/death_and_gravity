# Actual lapse insertion and local compact-support form

The original clock has a=(1+u^2)^2, h=(1+u^2)^3 and H=4u/(1+u^2).
Replay of S6.58 fixes the first mass variations, in the original
orthonormal frame, to

    Yt=m^2 alpha n, Ys=m^2 beta n,
    alpha=4/(9h), beta=28/(81h).

At the clock the unperturbed deviation Y is zero. Perturbations
of the frame or clock normal therefore do not add a first-order
off-diagonal Y. Changes of the scalar profiles' metric are also
irrelevant to this quadratic-in-Y term. The actual chain rule is

    D^j(alpha n)=sum_{r=0}^j binomial(j,r) alpha^(j-r) n^(r),

and similarly for beta. All H jets are replaced by derivatives
of the actual H, and spatial momentum follows k'=-Hk.
Freezing alpha,beta while differentiating would give a different
curved answer.

Let r=1+u^2, q=k^2, with k the dimensionless physical spatial
momentum. The full mass-insertion pole action is, modulo
compact-support boundary terms,

    integral a^3 [A_j n''^2+B_j n'^2+C_j n^2],

for covariant derivative orders j=0,2,4, in the common Euclidean
prefactor 1/(32 pi^2 epsilon). The exact coefficients are

    A0=B0=0,
    C0=-452m^4/(2187r^6),

    A2=0,
    B2=-832m^2/(6561r^6),
    C2=-16m^2[13q r^2+258u^2+78]/(2187r^8),

    A4=-128/(6561r^6),
    B4=-64[11q r^2+612u^2+90]/(19683r^8),
    C4=-16[101q^2 r^4+10q r^2(265u^2+283)
           +238230u^4+172980u^2-4770]/(98415r^10).

These are squared-frequency local residue coefficients, not a
renormalized finite inverse propagator.

For a raw quadratic density with coefficients c_ij multiplying
n^(i)n^(j), i<=j<=2, define D_w=D+3H. The integration-by-parts
reduction used by the executable derivation is

    A=c22,
    B=c11-c02-D_w(c12)/2,
    C=c00-D_w(c01)/2+D_w^2(c02)/2.

The Euler derivative of the raw-minus-reduced action vanishes
exactly, including Dk=-Hk. Each of the nine displayed coefficients
is independently compared with that reduction.

## Relationship to existing matching

The zero-derivative C0 already lies in the full S6.53 flat
potential; it must not be added again. Q2 and Q4 above are
quadratic mass-deviation terms missing from that checkpoint's
certified first-variation matching scope. The covariant definition
is a Taylor coefficient at fixed metric, not an instruction to
replace the existing pure-metric or linear-mass counterterms.

Every new quadratic operator has zero value and first variation
at Y=0. In particular it does not change the S6.60 selected
stress profiles or their first-variation cancellation. This
does not cancel any finite contact or retarded second variation.
The mixed metric/linear-Y terms, dimensional finite pieces,
initial covariance response and full integrated kernel must
still be included before a physical perturbation result.
