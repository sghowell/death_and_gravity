# Full conserved-preparation contraction

This proof treats the actual integrated-trace map in `actual_map.py`.
It retains the Einstein coefficient `1/(60 delta)`, the anomaly, the
named local coefficient, the complete nonlinear quantum response and
the prescribed conserved preparation source. The A.10 response estimate
alone is not substituted for this map. Smoothness, agreement with an
actual Hadamard state and propagation of the complete constraint require
the accompanying regularity and formulation arguments.

## 1. Named clock, ball, and provenance

Use the actual dimensionless conformal clock `x=eta/eta_star-1`, so
the last free slice is x=0. Write
`a_phys=A eta_star a`, `u=eta_star^2 U`, `h=a'/a`, and choose physical
epsilon=1, delta=10^-14. The prime in this note is d/dx. Thus
`kappa hbar=2880 pi^2 delta A^2 eta_star^4`. The ordinary radiation
normalization and the named length `lambda=2 sqrt(2) A eta_star^2`,
gamma=0, remain fixed.

Start at the old prepared label y0=5/2, and shift the future clock to
t=x-x0. The interval length is L=10^-10. The unknown is

    X=u'-ubar', W(t)=integral_0^t X(s) ds, u=ubar+W,
    X(0)=0, ||X||<=r=10^-6.

The entire history before x0 is unchanged. In particular the original
free-past plane-wave state is not reset at x0. The actual nonlinear
functional R is that of A.10 on this whole history.

A.8 supplies the following exact inputs through its source-hash-pinned
calibration. These are not new guessed state bounds:

    b0=20461/128, b1=61013499/8192,
    |ubar|<=delta b0, |ubar'|<=delta b1,
    |Sbar|<=delta S_Born,0+delta^2 R_0 <1.2*10^-10,
    |Sbar'|<=delta S_Born,1+delta^2 R_1 <4.2*10^-9.

Here the last two *exact*, unrounded values are recomputed in bounds.py
from the A.8 calibration, approximately 1.1499441689*10^-10 and
4.1620815964*10^-9. The decimal values in this sentence do not prove
the inequalities.

If cbar is the weighted actual old density defect, the A.8 bound gives

    |cbar|<=81 delta [14000+delta (8*10^14)+3/(2880*1922)]
            < C=12*10^-9.

The factor 81 is the coarse plateau bound a^4<=3^4. The pinned A.8
physical residual is in units eta_star^-4, not y^-4; no extra inverse
power of y is silently inserted here. The literal factor 1922 is the
one independently proved in A.8.

Let c(t)=chi_L(t)cbar(t), with 0<=chi_L<=1, chi_L identically one on
[0,L/4] and identically zero on [L/2,L]. Both |c|
and |c-cbar| are <=C. The known old cbar is smooth; it is a prescribed
function, not reevaluated on each unknown metric. No derivative bound
for chi or cbar is used in the C0 contraction.

The old plateau satisfies dy/dx=(1-delta/y^4)^(1/4)<=1, so its known
continuation stays in 2<y<3 for the whole future interval. The total
active history remains less than T=3. On the ball,

    |u| <= delta b0+Lr < B=2*10^-12,
    |ubar'|<V=10^-10,
    ||u'||_whole_history <= delta b1+r < M=10^-5.

The cutoff and geometry bound involve a closed finite interval only;
they make no claim about a much later cosmological evolution.

## 2. Geometric enclosure and pair estimates

At y0=5/2, the exact plateau formulas imply

    12/5 < a0 < 13/5, 3/8 < h0 < 5/12.

For example `a0^4=y0^4-delta` and
`h0=1/[y0(1-delta/y0^4)^(3/4)]`; raising the claimed positive
endpoint inequalities to fourth powers verifies them rationally.

Solve `h'=-u-h^2`, `a'=ha` with these unchanged data. As long as
1/3<=h<=1/2, `|h'|<=1/4+B<1/3`. A first-exit argument and
L/3<1/24 prove that h remains in this interval. Since a increases,
a>2. Also L<1/10 and

    a<= (13/5) exp(L/2) < (13/5)*(20/19) <3.

The exponential bound uses exp(z)<=1/(1-z) for 0<=z<1.
Thus every ball input produces a positive smooth-enough metric on
the whole interval; no a=0 or h=0 denominator is approached.

For two inputs let D=||X1-X2||. Their same-data Hubble difference is

    Delta h(t)=-integral_0^t exp[-integral_s^t(h1+h2)] Delta u(s) ds.

The exponential is <=1 because both Hubble functions are positive.
It follows pointwise that

    |Delta u(t)|<=tD, |Delta h(t)|<=t^2 D/2,
    |Delta log a(t)|<=t^3 D/6,
    |Delta a^2(t)|<=3t^3 D,
    |Delta a^-2(t)|<=t^3 D/12,
    |Delta d(t)|<=t^3 D/12, |Delta d'(t)|<=t^2 D/4.

The two scale inequalities follow from the mean-value theorem in log a
using 2<=a<=3. The local coefficient is the exact named one,

    d(a)=-19/60-log(a/2)/2, d'=-h/2, |d|<=17/30.

Indeed 0<=log(a/2)<1/2: exp(1/2)>1+1/2+1/8>3/2.

## 3. Exact auxiliary equation and removal of cutoff derivatives

Let q0=Sbar(x0)/a0^2 and q0'=(Sbar'-2h0 Sbar)/a0^2. The preceding
state bounds prove the deliberately coarse caps
`|qbar|, |qbar'|<=Qb=10^-8` on the old plateau. Put P=a^2 q'.
The complete forced trace integrates to

    P(t)=P0+integral_0^t [a^2u/(60delta)-u^2/4-h^2(u+h^2)/30] ds
         +8[c(t)/h(t)-c(0)/h0-integral_0^t c(1+u/h^2) ds],
    q(t)=q0+integral_0^t P/a^2 ds.

This is exact, since `(1/h)'=1+u/h^2`. It is not a bound for a
surrogate quantum source. Integration by parts includes *both* endpoint
terms. It avoids a false O(1/L) cutoff-derivative obstruction.

At X=0, the geometric terms coincide with the old solution, which has
the source cbar. Their source difference vanishes at t=0. Therefore

    E_P=8C[3+L(1+9B)],
    ||P-Pbar||<=E_P, ||q-qbar||<=L E_P/4,
    E_F=E_P(1+9L/4)<3*10^-7.

The last line bounds the actual fixed-point RHS at the center because
all its other difference terms vanish there.

For two arbitrary ball inputs the source c is identical. Pointwise,

    |Delta(1/h)|<=9t^2 D/2,
    |Delta(u/h^2)|<=9tD+27Bt^2D.

Termwise integration gives `||Delta P||<=K_P D` with

    K_E=3L^2/(40delta)+B L^4/(80delta),
    K_C=B L^2/4+L^2/240+(B+1/2)L^3/180,
    K_source=72C L^2+72CB L^3,
    K_P=K_E+K_C+K_source.

In particular the large Einstein scale is present in K_E; it was not
mistaken for a quantum-response constant. Numerically K_P is just above
7.5*10^-8, with that decimal used only for orientation.

Before using q/P caps, verify them without circularity:

    |P| <=9Qb+E_P+K_P r < Pmax=10^-6,
    |q| <=Qb+L Pmax/4 < Qmax=2*10^-8.

K_P does not depend on these proposed caps. Now mean-value integration
gives

    K_q=L K_P/4+Pmax L^4/48,
    K_G=K_P+9K_q+9Qmax L^2+3Qmax L^3,

where K_G bounds the change of `G=P+2ha^2q=a^2(q'+2hq)`.

## 4. Complete RHS Lipschitz constant

For the actual Wick functional

    S=R+4pi^2 D_EulerGamma[u]+d(a)u,

all common-past terms cancel. Derivative transfer has no old free-past
endpoint term, since the old potential was flat there. The equation for
X has the exact RHS

    F[X]=G-Gbar-(R'[u]-R'[ubar])-dX-(d-dbar)ubar'
                         -d'W-(d'-dbar')ubar.

A.10 supplies the nonlinear-response pair bound

    K_R=MT L^2[5/4+log(T/L)/2]+18M^2T^5 L exp(2MT^3).

Here T=3, M=10^-5 and z=MT^3=27/10^5. The rational bounds
`log(T/L)<32` and `exp(2z)<2` give the coefficient used in bounds.py:

    K_R_upper=MT L^2(5/4+16)+36M^2T^5L <10^-15.

For the logarithm, log(3)<2 and log(10)<3 follow from positive partial
exponential series; for the exponential use 2z<1/2 and exp(1/2)<2.
These are coarsenings of the actual same-history lemma, not an assumed
future state reset or a finite-frequency truncation.

Every local product-rule term is accounted for. In the order displayed
in F, their pair constants are

    dX:                   17/30 + r L^3/12,
    (d-dbar)ubar':         V L^3/12,
    d'W:                  L/4 + r L^3/4,
    (d'-dbar')ubar:        B L^2/4.

Thus

    K_local=17/30+L/4+B L^2/4+(r+V)L^3/12+r L^3/4,
    K_F=K_G+K_R_upper+K_local <567/1000.

All potential, geometry, auxiliary and prescribed-source changes enter
this estimate. The proof is on the closed ball in C0 with X(0)=0; no
derivative of X is required here. The response functional is the A.10
C1-valued nonlinear functional on such potential inputs. Agreement with
smooth actual states is addressed after the regularity argument, not
presupposed for every ball element.

## 5. Causal inverse and the finite ball

The required inverse is J_EulerGamma/(4pi^2), not J_EulerGamma itself.
From A.10's pole-plus-cut formula, EulerGamma>0, pi>3 and log(10)>2,

    k <= 2L/(1-L)+(2/9)sqrt(L)+4/log(1/L)
       < 2L/(1-L)+2/(9*10^5)+1/5 <21/100.

The first term retains the positive pole. The elementary inequality
pi>3 can also be checked without a transcendental approximation:
`1/(1+x^2) >= sum_(j=0)^7 (-1)^j x^(2j)` on [0,1], and four times
the polynomial integral is >3. The remainder is x^16/(1+x^2).
The bound log(10)>2 follows from e<3 and 3^2<10. The rational
certificate checks all displayed margins.

Consequently the *full* map has

    Lipschitz constant <= k K_F <3/25,
    center norm <= k E_F <63*10^-9,
    self-map norm <63*10^-9+(3/25)*10^-6
                   =183*10^-9 <10^-6.

This proves strict contraction and self-mapping on the stated finite
ball. It does not prove either condition for an enlarged interval; the
tests deliberately enlarge L and reject the corresponding gate.

The causal inverse maps C0 to C0 and sends every continuous RHS to a
function vanishing at the future origin. Because c=cbar on a fixed
initial neighborhood, the map also preserves that common zero
neighborhood along Picard iteration. This additional flat-start fact,
not a false generic C1-endomorphism assertion, is used by the separate
smoothness argument.

## 6. Verification and physical boundary

The source and state bounds are obtained from immutable A.8. The
nonlinear mode Lipschitz and normalized inverse bounds are the actual
A.10 lemmas. New arithmetic is exact rational arithmetic, with direct
Fraction tests of the integration constants and independent replay in
the parent verifier. The infinite-frequency, retarded-kernel and
Banach implications remain written mathematical proofs.

The preparation tensor is conserved and becomes zero before the future
half-interval. Its pressure depends on the unknown h; replacing it by
radiation pressure while switching c would violate conservation. The
construction does not alter the fixed ordinary radiation constant or
the quantum state's actual common-past integration constant.

This is a finite local preparation construction, not propagation of
the original A.8 geometry as an exact SEE solution. It does not transfer
the A.9 QSEI automatically, prove a cosmological singularity theorem,
provide long-interval shadowing, or close the realistic-field part of
P8(a). The old off-shell radiation past also precludes claiming a
globally compact-in-time preparation source on that entire past patch.
