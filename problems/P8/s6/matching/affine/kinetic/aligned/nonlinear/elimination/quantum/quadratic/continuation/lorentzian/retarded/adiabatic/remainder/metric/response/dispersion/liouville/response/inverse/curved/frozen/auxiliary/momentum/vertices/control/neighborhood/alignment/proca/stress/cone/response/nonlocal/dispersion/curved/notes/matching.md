# Actual curved scalar response and its fixed massive reference

This note proves the remainder for the new scalar current. It reuses
only the generic current-polynomial, finite WKB and time-primitive
arguments of S6.73, whose source and proof remain frozen upstream.
Its old lapse mass vertices, local matrix and full two-source inverse
are not used. The background mode equation and selected all-order
Cauchy state are unchanged by S6.81; the response vertices are new.

## Literal current, all dimensions and subtraction

The physical pure-scale current has the actual minimal-Proca
Hamiltonian form J_Z=(A_Z p^2+B_Z omega^2 f^2)/2. Its retarded bilocal
part is the current commutator. All canonical second-vertex,
normalization and fixed-profile contacts are retained. The latter
are local and do not change the following off-diagonal computation.

At high comoving momentum k, the actual current pair A-B, in physical
lapse/scale order and spatial dimension D, is

    b_T,D=(0,2(D-3)), b_L,D=(0,2(D-1)).

This is extracted from the full new minimal Hamiltonian before D=3,
not obtained by continuing only a physical scalar polarization.
The transverse multiplicity is D-1. Its highest pair square and its
first D derivative vanish at D=3. The longitudinal scale-pair
derivative is 2, and the derivative of its square is 16. These
evanescent contributions are kept.

At physical D=3, b_L,Z=4. With Q_w=64*pi^2 times the unscaled
3*delta p_Gamma, the leading retarded radial kernel is

    128 k^4 sin(2k Delta_sigma)/[a(t)^4 a(s)], sigma'=1/a.

The exact generic connected-current polynomial has k^4 coefficient
(A_t-B_t)(A_s-B_s) and k^3 coefficient
2i[b_s(A_t-B_t)A_s-b_t A_t(A_s-B_s)], where the b_t,b_s in this
last expression are canonical pumps, not the pair vectors.
They are kept as distinct functions. The diagonal single-scale
antisymmetric connection vanishes, but this does not set the whole
off-diagonal connection term to zero.

The Abel integral of k^4 sin(2k Delta_sigma) is
3/(4 Delta_sigma^5). Therefore the highest off-diagonal term is

    96/[a(t)^4 a(s) Delta_sigma^5].

For tau=t-s, the geometric ratio
tau^5/[a(t)^4 a(s) Delta_sigma^5] is smooth on the closed causal
triangle, equals 1 on the diagonal and has first lag derivative
-3H(s)/2. Four derivatives of 4/tau give 96/tau^5.
This matches the new scalar F_m high logarithm, not an old matrix
normalization.

## Finite selected-state expansion

For each desired finite derivative count, choose a sufficiently high
finite WKB order. The actual background pumps and potentials have
smooth rational dependence on k^2, with m>0 and a bounded away
from zero. The normalized variation-of-constants equation controls
the exact state against this reference by the finite Riccati
residual. The original all-order initial matching makes the error
arbitrarily inverse-power small at D=3, after the prescribed time
derivatives. It does not require replacing the selected state by
a reference vacuum or setting its initial interference to zero.

Expanding the finite reference amplitude and the residual phase
around exp(-ik Delta_sigma) gives a finite sum of smooth coefficients
times k^j sine/cosine terms, j<=4, and an integrable remainder.
Expand sufficiently far into negative powers for each desired time
derivative count. The finite low-k region is smooth and integrable;
the constrained canonical k=0 limit is the original one.

Only the undifferentiated common dimensional limit is required
near D=3. The unchanged W4(D) plus dimension-independent Borel
correction gives initial mixing O(k^-6). With the radial k^4
prefactor, its residual is O(k^(-2+|Re D-3|)), integrable for
|Re D-3|<1/4. At D=3 the stronger all-order matching supplies
arbitrary fixed time differentiability. No momentum derivatives
of the cutoff-summed state are needed.

A simultaneous diagonal time derivative of the phase inserts
k*(sigma'(t)-sigma'(s)); the bracket is O(t-s). Thus the accompanying
lag factor compensates the extra radial power. Deepening the finite
expansion for each fixed derivative count yields the same singular
class and an integrable differentiated remainder.

## No unspecified fourth-order local difference

The new S6.86 reference is the exact massive flat scalar F_m,
including its fixed finite coefficient -4 and its derived dispersion
constant. It is not replaced by a bare logarithm.

At frozen scale a0 the dimensional radial measure/amplitude has
a0^(-D-2) k^(D+1) dk. Under k=a0 p this is exactly p^(D+1) dp
for every nearby complex D. Hence the leading curved amplitude
agrees on the diagonal with the physical flat amplitude throughout
that dimensional neighborhood, including the first D jet.
There is no leftover log(a0) fourth-derivative contact.

For varying scale, their difference has a factor tau with an
analytic-in-D smooth coefficient. The causal leading power
tau^(-5+2epsilon) can carry delta'''' in its dimensional finite
part; multiplying by a diagonal-vanishing factor reduces the
worst power to tau^(-4+2epsilon), which carries at most delta'''.
This follows directly by subtracting the test function's Taylor
polynomial in the defining radial/lag integral: poles occur at
the Taylor powers, with no arbitrary local remainder. Matching
only at D=3 would not have excluded an evanescent delta''''.

The actual new order-four adiabatic coefficient is independently
extracted before freezing Hubble jets or taking D=3. All eight
two-current entries at every momentum equal
(A_i-B_i)(A_j-B_j)/32. Thus its highest input derivative is
precisely the new flat Taylor subtraction. The actual finite local
scale coefficient is independently -4 from S6.84; lapse fourth
coefficient is zero. The new massive reference retains all
finite low-momentum Taylor pieces, not just their high-k limit.

Second Hamiltonian contacts, output normalization and the fixed
profile are order zero. All other explicit local differences
have order at most three. Lower off-diagonal singularities are
causal finite parts tau^-n with n<=4, logarithms and an ordinary
integrable remainder. Their extension has just been fixed in
the same scheme as the actual response.

## Four primitives

After four zero-past primitives, tau^-n (n<=4) becomes a constant
multiple of tau^(4-n) log(tau), plus smooth causal polynomials.
For a smooth amplitude expand about the diagonal through n-1;
the remaining factor cancels the singularity. The high radial
k^-1 term has at most a logarithm, obtained by splitting at
k=1/Delta_sigma and integrating the oscillatory tail by parts.
Local derivatives through order three yield bounded polynomial
kernels after integration by parts. All prepared endpoint terms
vanish. The exact primitive identities and multiplier terms are
checked in code.

It follows for the actual pure-scale channel that

    I4 Q_w,scale = F_m(partial_u^2)+V_w,scale,

where every fixed diagonal derivative of V_w,scale is bounded by
C_j*(1+|log(t-s)|), with finite constants on the compact triangle.
No numerical values for these constants are claimed.

All quantum eta terms were proved local in coordinates.md.
Combined with the full actual tree inventory in tree.md, their
primitives give the asserted complete two-channel remainder.
The first row contains no logarithmic nonlocal term: after the
actual A eta'''' endpoint is separated, its local-primitive
kernel and first output derivative are bounded. This stronger
fact is retained for reconstructing the physical lapse.
