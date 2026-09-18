# One-block hard core for the three-soft-current subtraction

The executable current caps in this packet are premises proved by its
separate all-sector polynomial certificates, not external assumptions.
No original P8 closure is claimed.

## Object and exact scope

Use the unchanged original parameters, canonical action, positive Born
normalization A0=Am+AG, and real collision domain of S319: mass-one
external scalars, 5/4<=E<=2, W=a+b+c<=1/8, positive elastic t/u
transfers, fixed emitted directions and unit complex TT leaves.
Let H3 be the COMPLETE four-tree pure-soft temporal current of S317.
In the S319 maximal-soft partition, the one-block contribution has
all three emitted leaves in this one current. Its hard core has 47
trees; with the four internal soft trees the inventory is 188.

This is a term in the complete TEMPORAL tree reorganization. It is
not declared an independently gauge-invariant observable or identified
with188 isolated diagrams in a different gauge. The remaining cores
and the real-virtual/state matching are separate obligations.

Let R(Q;p) be the hard linear functional on an arbitrary symmetric
spatial input tensor, using the literal47-tree source with four
on-shell hard scalars and one possibly off-shell soft root Q.
No null, TT or multi-root Ward premise is imposed on that input.
Write the one-block term as R(H3). Define

 F=(abc/W)*H3,
 K=W*sqrt(kappa)*sqrt(rho)*R/A0.

Then its energy-scaled coefficient is G_block=K(F)/sqrt(kappa).
Both expressions use the same recoil map and analytic phase; W in
an analytic expression means the complex sum of energies, never an
absolute value. Real-center W and A0 are only bounding parameters.

## A global total-energy polydisc for the HARD core, not for H3

Fix a positive real center and take |z_i-w_i|<=eta*W for all three
energies, eta=10^-13. Then |delta Q0| and the spatial Euclidean
norm of delta Q are at most3*eta*W<10^-12*W. The recoil square-root
estimates in S313 notes/analytic.md allow perturbations at least this
large (there the two-energy total perturbation was2*10^-12*W).
The S320 extension of those estimates depends on the total perturbation,
not the number of leaves. Square-root branches, massive shells and the
analytic phase therefore continue with the same strict margins.

The hard core contains only ONE soft block. Any external scalar-chain
cut containing radiation therefore contains ALL of Q, not a small
proper subset of a,b,c. Its real gap is at least3W/8 and its complex
perturbation is below100*10^-12*W, leaving inverse below4/W.
Hard mixed cuts retain magnitude>(tau+W^2)/600000; timelike and heavy
cuts retain45/16 and n/2. These are the S313/S320 bounds on total-Q
assignments and hold even when the one block is timelike rather than
null. The S319 core proof explicitly allows composite off-shell blocks.

Crucially R has no pure-soft root or child invariant in a denominator.
An off-shell root is amputated. With one block there cannot be another
nontrivial pure-soft branch inside the hard core. The Q^2=0 obstruction
of S320 concerns H3, which has been factored OUT here; no global-disc
claim is made for H3 or the full three-real amplitude. All remaining
hard-core denominators are nonzero on a neighborhood of the closed
polydisc. K is therefore operator-valued holomorphic there.

## Explicit operator envelope from the existing core majorant

The S320 complex hard-core functions, before composition with soft
currents, are

 L=(1-1024x)/(1-5120x), D=1/(1-4x),
 H=(1-1024x)/(1-3072x), E=1024/(1-1024x)^2,
 V=120000000*1024*(2/(1-32x)^3-2),
 Fm=L^4*(D+(3/2)*D^2*H),
 FG=(3*600000/8)*L^4*E^2/(1-V).

Their linear coefficients are
 cm=44048,
 cg=5566277620447838208000000.
Direct symbolic differentiation and independent coefficient convolution
agree. The vertex norms in S319 notes/vertices.md apply to arbitrary
complex tensors, so these are multilinear operator estimates, not only
bounds at physical TT test inputs. The one-block charge4/W is retained;
an uncharged block is enlarged by1/W as in the original proof.
Thus the unscaled core has norm at most
 (n^2*cm+cg)/(sqrt(kappa)*W) relative to A0.
Since |sum z_i|<2W and |sqrt(rho)|<2, the normalized operator obeys

 norm(K)<=B0=4*(n^2*cm+cg)<10^8*n^2+10^28.

Cauchy in the three energy variables now gives, for every square-free
multiindex alpha,
 norm(partial_alpha K)<=B0/(eta*W)^|alpha|.
The factorial is1 for these multiindices. No physical matching constant
has been selected; B0 and eta are explicit conservative bounds derived
from the fixed source and its existing complex-gap inequalities.

At any energy face with W>0 the same hard denominators remain nonzero,
so K extends continuously there. A common limit at W=0 is unnecessary:
K is bounded while F=O(W^2), so K(F) has the zero limit at the origin.

## Product and integrability statement

The replay must prove the all-angle caps C1,C2,C3<=10^6
in the hierarchy
 norm(F_i)<=C1*W/kappa,
 norm(F_ij)<=C2*W*S/kappa (i not j),
 norm(F_abc)<=C3*S/kappa,
 S=1/(a+b)+1/(a+c)+1/(b+c),
and compatible faces. S318 already supplies
 C0=2*(2*10^10)^2*3!=4800000000000000000000
in norm(F)<=C0*W^2/kappa for unit leaves.
The exact eight-term product rule and2*W*S>=9 give

 norm(partial_abc G_block)<=C_block*S,
 C_block=B0/kappa^(3/2)*
 [C3+3*C2/eta+(2/3)*C1/eta^2+(2/9)*C0/eta^3].

At the unchanged n=10^200/512+2 and kappa=10^800, exact Fraction
arithmetic proves C_block<10^-738 using the displayed executable caps.
The exact polynomial certificates, not this arithmetic alone, prove
the current caps for all polarizations in both independent sectors.
The compatible-face fundamental theorem then gives
norm(Delta_a Delta_b Delta_c G_block)<=10^-738*J(a,b,c), where
J=c*I(a,b)+b*I(a,c)+a*I(b,c) and
I(a,b)=(a+b)log(a+b)-a log(a)-b log(b).
The known integrals of J/(abc) and J^2/(abc) then apply to this class.
They do not establish the other4,928 temporal-tree terms, any physical inclusive
probability, an all-N sum, a loop remainder, Regge control or P8 closure.
Here5116-188=4928 is only the complementary inventory, not a claim
that a single common subtraction bound already controls it.

## Independent original-parameter calibration

Two real recoil configurations and their exact Gaussian-rational
continuations use unchanged original n,g,C,kappa. Real points agree
with S300 recoil and retain positive Born pieces. Complex points are
inside the existing relative tube; this calibrates the source, not
the larger global hard-core tube proved above. At all four points,
the47 hard-root source is conserved, all six spatial root orientations
agree with independently rooted scalar amplitudes, and the Hermitian
operator norm satisfies B0. The complete four-tree current contracted
with that root agrees with the collapsed47 core and separately filtered
188-term complete temporal class. No omitted source is relabeled full.
