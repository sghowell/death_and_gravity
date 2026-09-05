# Coupled M1 finite-time tree estimate

The theorem is exactly the [operational contract](../FORMULATION.md).
All constants are deliberately conservative sufficient bounds. The
recurrence and final exact fractions are in `majorant.py` and the replay
certificate. The proof below explains why those arithmetic bounds cover
the claimed continuum, time and field domain.

## 1. Fixed local units and the complete free-mode band

Set L=10^11, U=10^12, s=ell(t)/ell0, and a0=1. On
|t-t0|<=ell0/100, S5.7 gives

    s in [99/100,101/100],
    |x-x0|<=1/99, |log(a/a0)|<=4/99,
    95/99<=a/a0<=99/95,
    E(t)/E(t') in [7/11,11/7].

We may loosen s and a to [1/2,2]. The prescribed one-chart selector stays
admissible throughout the interval. Every external or internal tree mode
has physical momentum, in fixed ell0 units, between L/2 and 4U. Thus its
local q obeys

    q>=L^2/16>10^20, q<=(8U)^2<qmax=10^27.

Every degree-at-most-four Fourier-mask physical derivative is bounded by
K=10^13. Required nonzero York inverses have k^2>=(L/2)^2. The full
four-leg total mask has zero momentum, but the trace/TT argument in S5.6
proves that its fourth-order York solution cannot enter H4: only the first
three perturbative constraint orders are used. Homogeneous dynamics are
not thereby solved.

For the positive matrix R=W(left)^(1/2), the initial scalar mode matrices
U0=(2R)^(-1/2), P0=-i*(R/2)^(1/2) satisfy

    U0*U0^dagger-conjugate(U0)*U0^T=0,
    P0*P0^dagger-conjugate(P0)*P0^T=0,
    U0*P0^dagger-conjugate(U0)*P0^T=i*I,
    (P0^dagger*P0+U0^dagger*W(left)*U0)/2=R/2.

Here conjugate denotes entrywise conjugation. These are
canonical CCR and an initial covariance choice, not the assertion that
Hnorm=E-P^T*Omega*Y is diagonal or minimized by this state. In particular
Ydot0=P0-Omega(left)*U0. The real evolution matrix

    Aham=[[-Omega,I],[-W,-Omega]]

satisfies Aham*Jcan+Jcan*Aham^T=0 for the canonical antisymmetric Jcan,
so exact evolution preserves the CCR. `free.py` independently verifies
these identities for a generic coupled matrix and a rotated positive
initial spectrum. No time-dependent instantaneous eigenbasis is used.

The band and S5.7 potential bound give lambda_max(R)<5U and W>=L^2*I/8
in fixed units. Each initial scalar column therefore has E_j<5U/2 and,
throughout the window, E_j<(11/7)*(5U/2)<4U. It follows that

    |Y_column|^2<64U/L^2<1,
    |P_column|^2<8U<(10^14)^2.

The tensor initial data obey the same bounds. Conjugate columns obey them
as well. Each end of an internal contraction is one of these exact mode
columns. Two scalar and two tensor columns are required; a scalar
off-diagonal propagator is not deleted by calling the principal speeds
equal. The independent Taylor enclosure for exp(4/11) checks the slack in
the rational energy-growth estimate.

## 2. Exact canonical seeds, including the momentum boundary

Write alpha_bar=Tbar^T*Tbar. S5.7 supplies

    (1/100)I<alpha_bar<512I,
    ||Tbar||<23, ||Tbar^-1||<10.

With Lw=diag(s^-delta,1), delta=1 in gamma and zero in unitary, the
fixed-unit T is Tbar*Lw. Hence ||T||<46, ||T^-1||<20; also
a^-3/2<3. The exact inverse transformation gives

    Q=a^-3/2*T^-1*Y,
    p_chart=a^-3/2*T^T*(P-Sfixed*Y).

The symmetric momentum boundary cannot be dropped. Applying S5.7's exact
denominator/coefficient method to **z*Sbar**, rather than Sbar, gives a
uniform CS with ||Sbar||<=q*CS. The two exact chart bounds appear in the
certificate; their common upper bound is

    CS=42828424431/2000000.

Since Sfixed=Sbar/s, ||Sfixed||<=2*qmax*CS. With Wmode=10^14, each
column obeys

    |Q|<60,
    |p_chart|<138*(Wmode+2*qmax*CS).

The gamma inverse is zeta=p_gamma/(2*qfixed), p_v=-2*qfixed*b, where
qfixed=q/s^2 lies between L^2/4 and (4U)^2<qmax by the physical momentum
band, a stronger upper bound than a generic s-only conversion. The same gamma condition used
to normalize the free system holds at all external and internal waves.
In both charts the required mixed inverse is

    old_p=p_v+3*l_fixed*matter,
    pi_chi=l_fixed+P_m+3*l_fixed*zeta.

Here |l_fixed|<=1/5. Both shifts are retained. They are not substitutes
for the kinetic off-diagonal entries. Normalize each real TT basis tensor
to E:E=1. Its entries have modulus at most one; the tensor coordinate and
momentum follow the volume-normalized free map, including the 3H/2
generator, and satisfy weaker bounds than the scalar ones above.

Exact rational comparisons show that the common **terminal phase seed
B=10^34** bounds every normalized-mode contribution to zeta, gamma_ij,
old_p, PiTT_ij, matter, and delta(pi_chi). No separate independent seed
for only an unmixed matter mode is assumed. The new time generators of
the linear canonical map are quadratic; H3 and H4 transform by this
linear phase-space substitution. Free-velocity substitution into the old
unnormalized velocity kernels would miss the q-dependent S terms.

## 3. Positive coefficient induction through degree four

Use a nonnegative series in a formal degree marker e, truncated at degree
four. At each degree its coefficient bounds the sum of absolute values
of all labelled Fourier-mask coefficients of that degree. Convolution
therefore bounds products. A derivative multiplies the bound by K, and
discarding only degrees known to vanish is valid. Four labels per terminal
field, b=4B*e, bound both the cubic and quartic constructions.

Let r_g=3b bound each entry of 2*zeta*I+gamma. The following positive
series dominate the full metric geometry:

    metric <= 1+r_g,
    inverse <= 1+sum_(n=1)^4 3^(n-1)*r_g^n,
    det_remainder <= [6*(1+r_g)^3]_(degree>=1),
    volume <= absolute_binomial_(1/2)(det_remainder),
    inverse_volume <= absolute_binomial_(-1/2)(det_remainder),
    connection <= (9/2)*K*inverse*r_g,
    curvature <= 9*inverse*(6K*connection+18*connection^2).

The determinant has exact background one, not six. The factor six only
bounds its positive-degree permutation terms. Absolute binomial
coefficients through degree four give coefficient bounds on the formal
Taylor jets; they do not assume convergence for arbitrarily large
quantum-field values. Index counts cover all three spatial directions.

The free metric momentum density before York corrections is bounded by

    Pi_free <= 8+18b,
    Dchi <= 1/5+b.

Indeed the exact terms are -H*(1+zeta)*I+old_p*I/6+PiTT+H*gamma with
|H_fixed|<=8, while Dchi denotes the full matter density. For each
n=1,2,3, take the degree-n part of

    source <= 3K*Pi+9*connection*Pi
              +(3/2)*K*Dchi*inverse*b.

The final term is the **matter source** in all three spatial constraints;
it is absent in D-only M0 and cannot be omitted here. The exact York
inverse is (I-kk^T/(4k^2))/k^2. Each vector component is bounded by
2*source/kmin^2, and each York tensor component by 4K times that. Add

    correction_n=(8K/kmin^2)*source_n, kmin=L/2,

at degree n to Pi before forming the next source. This majorizes the
actual triangular nonlinear constraint recursion, not an inverse of a
sampled matrix. The certified nonexceptional domain supplies every
denominator lower bound.

Set mixed<=3*Pi*metric*inverse_volume. The five invariant bounds are

    sigma <= [3*mixed]_(degree>=1),
    rho <= curvature,
    eta <= [Dchi*inverse_volume]_(degree>=1),
    shear2 <= [12*mixed^2]_(degree>=2),
    z_matter <= 9K^2*inverse*b^2.

The removed backgrounds and linear shear/gradient terms vanish by the
exact S5.6 identities. In particular eta is the perturbation of the
**density divided by volume**, not just the free matter momentum, and
z_matter uses the full inverse metric.

For a compact invariant monomial with powers (p,r,e,s,z) in these five
invariants, the fixed-unit coefficient carries

    (ell/ell0)^(p+2r+e+2s+2z-2).

The exponent lies in [-2,6] through the retained order, so its modulus
is at most 64. Insert all recorded M1 stationary Hamiltonian coefficients
with their S5.5 uniform absolute bounds multiplied by this factor. The
action's overall a^3 is at most 8. The -2H*Pi:metric boundary needs
144*Pi*metric inside that factor, because M1 has |H_fixed|<=8. The
explicit Adot, gamma and clock-shift terms are at most quadratic in
this linear phase parametrization; a harmless positive majorant can
retain them without using a cancellation.

Call the resulting positive Hamiltonian series h. The finite recurrence
gives exact rational bounds on its coefficients. Set

    B3=720*h_3, B4=720*h_4.

The redundant 6! factors safely include normalized Bose occupation and
tree Wick-label combinatorics up to four external legs and two ends of
one internal contraction. Each individual labelled kernel is already
bounded by the sum-norm coefficient; repeated fields have distinct labels.
No cancellation of a frozen vertex or asymptotic principal symbol is used.

The prior construction is a finite rational tensor expression in momenta
and polarizations away from the declared poles. Its identities and this
absolute induction extend from exact rational fixtures to real momenta
and normalized real TT polarizations by continuity, or directly by the
same index inequalities. Rational-angle sampling is not the proof of
coverage. A complete expanded symbolic quartic kernel is not needed for
this norm theorem.

## 4. Tree time ordering and the continuum operator bound

S5.5's verified scale restoration gives an overall (M*ell0)^2 in fixed
patch coordinates. Normalizing each perturbation canonically therefore
gives H3/(M*ell0) and H4/(M*ell0)^2. The additional linear canonical map
is independent of this overall factor; every time generator was already
retained in the quadratic/free analysis. The matter amplitude scales with
M as required by the fixed backreacted witness, not as a fixed external
spectator when tau is changed.

Evaluate the interaction Hamiltonians on the exact free columns. The
window's fixed-unit duration is 1/50<1. The cubic entry is bounded by
B3/(M*ell0). At quartic order there are three external 2+2 partitions,
four internal columns (two scalar and two TT), and two time orderings.
The loose coefficient 24 therefore gives the entrywise tree bound

    [B4+24*B3^2]/(M*ell0)^2.

Both ends of a coupled internal scalar contraction can contain both
original species; those contributions are already included in the
terminal seeds for each free column. An extra diagonal-propagator
assumption would be invalid. Momentum conservation fixes the one internal
tree wavevector. No loop momentum integral is hidden in the factor 24.

The phase Hamiltonian route includes H4 and H3 exchange consistently.
There is no need to add a separate velocity Legendre contact or a contact
obtained by differentiating an already time-ordered propagator. Those
representations require their own matching cancellations. Quantum
ordering differences, contractions within one vertex and other loops are
outside the stated tree order, not assumed numerically small.

Spatial homogeneity produces the total-momentum delta. First decompose
the continuum one/two-particle spaces into fibers of P; the map
(k1,k2)->(P,k1) has Jacobian one. The one-particle fiber has four labels.
The two-particle fiber has a bounded relative-momentum domain and at most
sixteen ordered species pairs. If

    N=4*(4U+1)^3,

their respective measures are bounded by 4 and 4N. N is an intentionally
loose **measure** bound for d^3k/(2*pi)^3, not a physical state count.
The common D=8N^3 exceeds both measures and their geometric mean.
For an entrywise bound b0, row and column integration or direct
Cauchy--Schwarz gives operator norm at most b0*sqrt(mu_in*mu_out)<=D*b0.
Normalized Bose symmetrization changes no continuum measure and is
covered by the finite combinatorial allowance above. The bound is uniform
in P, hence bounds the direct-integral operator on normalized wavepackets.

Apply the hard-transfer mask to the quartic-order kernel before this
estimate. Removing entries can only decrease this absolute majorant;
the omitted forward/soft entries themselves have not been bounded.
Vacuum production is excluded, so no attempt is made to treat a
total-momentum delta as an ordinary function or to take a finite physical
box and discard its homogeneous gravitational constraints.

Set C3=D*B3 and C4=D*(B4+24B3^2). Since ell0>=tau, a common sufficient
choice is

    M*tau >= 1000*C3,
    (M*tau)^2 >= 1000*C4.

The replay records C3, C4, and the first decimal power meeting both exact
inequalities: **M*tau=10^324**. This finite, very large bound closes the operational gate
at every center time. It is not a necessary duration or a useful estimate
of where perturbation theory actually fails.

## 5. Evidence and remaining gates

The input physical Hamiltonian and exact coupled normalization remain
pinned and immutable. Tests verify the independent positive-series
convolution/binomial algebra, the canonical initialization/CCR, every
scale and Schur inequality, and controls that omit matter terms or shrink
the domain below the established free band. Exact output and scope
mutations are rejected by read-only replay. This is a certificate plus
written induction/ODE/operator proof, not Lean formalization or a claim
of independent external review.

The result supplies the M1 counterpart of the scoped D-only hard tree
criterion. It does not bound higher trees, loops or perturbative remainders,
compute radiative protection of the exceptional principal relation,
establish a global vacuum or asymptotic S matrix, control infrared or
nonlinear backreaction, or match the tube to an accepted UV parent. These
remain genuine additional research questions; full P8 is not completed.
