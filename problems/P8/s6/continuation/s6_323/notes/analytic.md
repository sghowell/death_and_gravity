## Uniform central sources without fake soft poles

Remove the external-chain factor L(x)^4 from S320's nonnegative
complex hard-core majorants; do NOT remove central soft attachments.
There is no scalar soft propagator to charge, and no pure-soft
current denominator. With
 D=1/(1-4x), H=(1-1024x)/(1-3072x),
 E=1024/(1-1024x)^2,
 V=120000000*1024*[2/(1-32x)^3-2],
the central coefficients are
 m_r=r![x^r](D+3D^2H/2),
 g_r=r![x^r]((3*600000/8)*E^2/(1-V)),
 C_r=n^2*m_r+g_r.

Their exact first four values are
 m=(5/2,3088,18923696,174400047744),
 g=(235929600000,5566277616582367641600000,
 262649930314558756486872922324992000000,
 18590067899792945264953972626291745126348750848000000).

Thus |H_R|/A0<=C_r/kappa^(r/2) on the hard-only global energy tube.
For the uses below only one energy may vary, or all three vary by
at most eta*W, eta=1e-13. Total momentum perturbation is at most
3eta*W, strictly inside S313's recoil/hard-cut estimates.
There is no relative-energy or collinearity premise for this central
operator: all remaining denominators are hard mixed, timelike or heavy.
The first hard inverse has magnitude<600000/(delta+W^2); each
further EH vertex plus hard inverse has the already derived120000000
paired momentum budget. Scalar endpoint momenta have components<3.
These are the same vertex/path arguments as the frozen majorant,
with the external chains absent, not an inference from a sampled norm.

## Five contributions needing no new forward Taylor theorem

Keep the analytic phase |sqrt(rho)|<2 and normalized singleton
|S_i|<32. Exact face limits are taken at fixed elastic parameters.

Connected triple on one line:
its all-angle bound from the generic13-tree proof is
 |U_123|<1e24*abc*S,
 S=1/(a+b)+1/(a+c)+1/(b+c).
There are four lines. Multiplying by H_empty and the phase gives
coefficient8*C0*1e24/kappa^(3/2) timesabc*S. Every face is zero.

Pair remainder plus singleton:
for a fixed pair there are16 line placements. Its total modulus is
below1024*C0*1e11*ab/kappa^(3/2). It vanishes on the a/b faces.
Continue c alone on radius eta*W. Cauchy gives c-derivative below
that coefficient divided by eta*W. Integrating c from0 and using
1/(a+b+t)<=1/(a+b) gives a rectangle bounded by that coefficient
timesabc/[eta*(a+b)]. Summing the three pairs givesabc*S, not3S.

Two central leaves:
for each central pair there are four possible singleton lines.
The same reasoning gives coefficient256*C2/(eta*kappa^(3/2))
timesabc*S, since the prefactor is exactlyab and all other factors
are hard-analytic in c after cancellation of the singleton pole.

Three central leaves:
all faces vanish because ofabc. The bound is2*C3*abc/kappa^(3/2).
SinceS>=1 on W<=1/8, enlarge it to the same coefficient timesabc*S.

One central leaf with connected external pair:
there are four pair lines. The full t_jk bound gives coefficient
8*C1*2e10/kappa^(3/2) timesabc*S, already zero on all three faces.

All useabc*S<=J, since integrating the decreasing kernel
1/(x+y) over a rectangle gives at least its endpoint value.
None of these estimates requires a full current global-W disc.

## Pair leading term plus singleton

Let ell_ij=L_ij/u. It is odd in p and, after division by
M=1e14*ab/u^2, its magnitude and massive-spatial gradient are
strictly below32 and2048. These are within the normalized-current
bounds in the S322 double-external forward-grouping proof.

The actual expression is the sum over4*4 pair/singleton line choices,
each with its exactly shifted Born hard kernel. It has NO ordering
correction, because it is defined by the exact scalar-line cumulant
decomposition. Its assigned-stress, hard-denominator, recoil and
phase comparison are a subset of the explicit S322 budgets.
Grouping the two legs of each mixed channel uses oddness of ell,
not a transverse composite-graviton assumption.

After subtracting the Born product, which is independent of c, the
remainder has modulus below
 1e14*B2*(ab/u)*W/kappa^(3/2),
 B2=1e40*n^2+1e60.
Only c needs continuation, with radius eta*W; a,b are fixed.
Its c derivative is below1e14*B2*ab/(eta*u*kappa^(3/2)).
The a/b faces vanish, so integration in c gives
1e14*B2*abc/[eta*u*kappa^(3/2)].
The arithmetic probe uses the still larger2e14*B2/eta^3 coefficient.
This leaves the two grouped leading sectors treated below, using
the endpoint identities in notes/grouping.md.

## A separate analytic neighborhood at zero radiation

For a mixed channel take sigma=1e-14 and the origin polydisc
|z_i|<=rho=sigma*t. Born t<4, so this radius is tiny even at the
largest angle. No positive-energy assumption is needed inside
this origin polydisc; gaps are proved directly from Born data.

Let Z0=sum z_i, v=sum z_i*n_i. Then |Z0|,||v||<=3rho and
h=E-Z0/2,
Eprime^2=h^2-v.v/4,
rprime^2=Eprime^2-1.
The difference of Eprime^2 from E^2 is at most6rho+4.5rho^2<7rho.
Using E>=5/4 and r0>=3/4, the original square-root branches have
|Eprime-E|<7rho and |rprime-r0|<14rho. The recoil formula
 p0=h-/+rprime*(u.v)/(2Eprime),
 p_space=+/-rprime*u+
 [+/-(rprime/Eprime)*(u.v)/(4(h+Eprime))-1/2]*v
therefore gives |dp0|<20rho,||dp_space||<40rho.

Every hard momentum changes componentwise by<100rho. For a mixed
Born denominator D0=-t^2,
 |D-D0|<800rho*t+40000rho^2<801sigma*t^2<t^2/2.
Both hard momenta and a possible central soft momentum have
Euclidean norm<4t. Thus an EH3 step has momentum-square/next-gap
ratio<32, already covered by the120000000 envelope.
Timelike and heavy gaps stay strictly nonzero as well.

Every Doppler factor stays>1/8. Each paired current is bounded by
1024*(t+40rho)<2e4*t. Repeating the exact endpoint grouping gives
the origin-disc estimate
 |F_m(z)|<=C*t^m
for each mixed channel. For the compact contact/heavy/timelike
families the same estimate holds with t=1.
The phase remains analytic and below2, already included in C.

This second disc is used only for these hard-only grouped kernels.
It neither asserts a full-current origin disc nor contradicts
the proper soft-cluster poles recorded in S320.

## Two-regime Taylor estimates, with uniform angular constants

For the physical-center disc, sum|z_i|<2W.
First suppose W<=sigma*t/4. Then z is inside the half-sized origin
disc. Cauchy about any point on the segment0->z gives:
- for m3, every second derivative is bounded by8*C*t/sigma^2;
- for m2, every first derivative is bounded by2*C*t/sigma.
Taylor's integral remainder and t<4 give, conservatively,
64*C*W^2/sigma^2 and16*C*W/sigma, respectively.

If W>sigma*t/4, then t<4W/sigma. The physical-disc estimate and
(t+W)^2<=2(t^2+W^2),delta<=t^2 give
 |F3(z)|<=2C*t^2*(t+W)<=160C*W^3/sigma^3,
 |F2(z)|<=2C*t^2<=32C*W^2/sigma^2.
The origin bounds give
 |F3(0)|<=64C*W^3/sigma^3,
 |sum z_i*d_i F3(0)|<=32C*W^3/sigma^3,
 |F2(0)|<=16C*W^2/sigma^2.
Since W<=1/8, these imply enlarged remainders
256C*W^2/sigma^3 and64C*W/sigma^2.

Sum the two mixed, timelike and matter families. The choices
 B3=1e4*C/sigma^3,
 B1=1e4*C/sigma^2
strictly dominate both regimes and all four families.
The comparison is made on the entire closed physical-center disc,
not only along an equal-rate real soft limit.

## Anchored rectangles

The constant and linear origin polynomial ofT3 is annihilated by
the threefold anchored rectangle. Its remainder is holomorphic
and bounded byB3*W^2 on radii eta*W. Three-variable Cauchy gives
 |partial_abc T3|<=B3/(eta^3*W).
Use1/W<=S and the compatible closed soft faces to obtain
 |rectangle(T3/kappa^(3/2))|<=B3*J/(eta^3*kappa^(3/2)).

For w_i*T1_i, the w_i face is zero. At fixed positive w_i,
the two-variable Cauchy bound from T1_i-T1_i(0) is
 |partial_j partial_k T1_i|<=B1/(eta^2*W).
Its twofold integral is bounded byI(w_j,w_k), because
1/(w_i+y+z)<=1/(y+z). Summing the three central choices gives
 B1*J/(eta^2*kappa^(3/2)), not three times that coefficient.

The other six contributions are bounded in part I. All define
compatible faces: hard-only kernels extend ordinarily at fixed
positive Born transfers, the line connected functions have uniform
zero faces, and the pair-leading function has its ab/u vanishing.
Trimmed rectangles and local integrability justify the zero-face
limits without treating an original propagator pole as a point value.

The exact original-parameter arithmetic places the sum of all eight
3767 coefficients, plus the retained nonsingleton coefficient, below
1e-670. This gives the full5116-tree AMPLITUDE rectangle bound
 |rectangle(abc*sqrt(rho)*M_3/A0)| <1e-670*J.

The amplitude argument alone does not establish an inclusive
probability. Its separate defined signed-measure consequence is proved
in notes/measure.md. This still does not establish a real-virtual
match, all-N sum, interacting state, loop errors,
absolute complex Regge, common-parent bounce or original P8 closure.

## Scope of the reused pair-leading comparison

Only the double-external leading part of S322 is used for L_ij:
there is no claim that this cumulant is an arbitrary off-shell tensor.
At fixed a,b, ell=L_ij/(a+b) divided by1e14*ab/(a+b)^2 has modulus
below32, spatial gradient below2048 and exact oddness on both massive
branches. These are precisely the current hypotheses of the owner
grouping and stress/denominator comparison in that part of S322.

The assigned momentum Q_pair=a*n_i+b*n_j is independent of the
endpoint owner. The Born comparison is independent of c. The S322
regular, ordering and trace-correction budgets are nonnegative
overcounts and need not be used, so its full B2 remains a valid larger
bound. Only c varies for the Cauchy derivative; a,b and its positive
normalization are held fixed. The resulting ab/(a+b) factor supplies
the a/b face limits. No a/b derivative of this normalization is taken.

All triangle estimates are for complex moduli on the relevant closed
disc, with strict margins. The output retains conservative exact
integer/rational budgets; the displayed topology, factorization,
branch and two-regime arguments supply their analytic hypotheses.
