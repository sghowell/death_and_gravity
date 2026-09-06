# Actual nonlinear Wick response in integrated-trace regularity

This note proves an actual two-potential response estimate. It does not
establish a semiclassical fixed point, a compatible initial-data set, or a
smooth/Hadamard solution of the SEE. Earlier A.1--A.9 files are unchanged.

## 1. Domain, normalization, and what is being extended

Use actual conformal time, translated so that a free plane-wave starting
slice is at 0. This may be later than the originally specified Cauchy
slice, provided the potential vanishes between them: free propagation
then changes only the irrelevant overall mode phase. Let `T>0`, and
let `U,V` be real functions in
`C1[0,T]`, with `U(0)=V(0)=0`. They need not be nonnegative. Assume

    ||U'||_infinity, ||V'||_infinity <= M,
    D = ||U'-V'||_infinity, B = M*T, z = M*T^3.

Then `||U||,||V||<=B` and `||U-V||<=T*D`. The quantum field is still
massless and minimally coupled. Its normalized exact modes solve

    w_k'' + (k^2+U)*w_k = 0,
    w_k(0)=1, w_k'(0)=-i*k, v_k=w_k/sqrt(2*k),
    F_k=|w_k|^2=2*k*|v_k|^2.

The value at `k=0` is the regular Volterra limit for `w_k`; the factor
`1/sqrt(2*k)` is kept only when forming the physical field. The radial
measure `k dk` below makes this massless infrared behavior integrable.

The functional being bounded is the nonlinear part

    R[U](t) = integral_0^infinity k*(F_k[U]-1-F1_k[U]) dk.

For smooth potentials flat in a past neighborhood, this is the nonlinear
term of the actual transported-state Wick square used in A.7:

    a^2 W = -hbar*K_{a,lambda}[U]/(8*pi^2) + hbar*R[U]/(4*pi^2).

Neither `R` nor the following comparison changes the named subtraction
length or finite curvature coefficient. The local terms and the
singular *linear* retarded logarithmic term remain outside `R`.

For general `C1` potentials, the Abel prescription and convergent series
below explicitly define a `C1[0,T]`-valued nonlinear functional. It agrees
with the actual Hadamard Wick contribution on the smooth, past-flat
subset. No uniqueness of extension is inferred from that subset:
smooth potentials flat near the initial slice are not dense in the full
`C1` domain when `U'(0)!=0`. Uniqueness here refers only to the value of
the specified Abel/series definition. Nor does this definition assign a
smooth Hadamard spacetime to every element of a `C1` ball. A metric
solving `a''+U*a=0` for such a potential has only `C3` regularity a priori.

## 2. Exact Dyson expansion and endpoint terms

Set `w_0(t)=exp(-i*k*t)` and

    w_n(t) = -integral_0^t sin(k*(t-s))/k * U(s)*w_{n-1}(s) ds.

These are the terms of the actual convergent Volterra/Dyson solution,
not modes of an auxiliary quantum state. Define the degree-n polynomial

    F_n[U] = sum_{j=0}^n w_j[U] * conjugate(w_{n-j}[U]).

In particular `F_0=1` and

    F_1[U](t) = -(1/k)*integral_0^t U(s)*sin(2*k*(t-s)) ds.

Write `S_k q(t)=integral_0^t sin(k*(t-s))/k*q(s) ds`. Integration by
parts, including the lower endpoint, gives

    (S_k q)' = S_k(q') + sin(k*t)/k*q(0).

Consequently the boundary in differentiating `w_n` is

    -sin(k*t)/k * U(0)*w_{n-1}(0).

For `n=1` this is `-sin(k*t)*U(0)/k`; for `n>1` it vanishes because
`w_{n-1}(0)=0`. Induction, with the explicit hypothesis `U(0)=0`, gives

    w_n' = -i*k*w_n + D_U w_n[U](U'),
    F_n' = D_U F_n[U](U').

Here `D_U` is the derivative of the n-homogeneous multilinear polynomial:
it replaces one potential factor by the supplied direction and sums all
insertion positions. The free `-i*k` and `+i*k` phases cancel in every
product contributing to `F_n`. No further integration by parts and no
second or higher derivative of the potential is used. The inserted
direction `U'` need not vanish at the original endpoint.

For the full smooth past-prepared geometric interpretation, all jets at
the original join vanish. That stronger hypothesis is not needed for
this one-derivative nonlinear estimate. Omitting the displayed endpoint
for a nonzero initial potential is an explicit negative control.

## 3. The complete second-order contribution

Direct multiplication of the first three mode terms yields

    F_2(t) = (1/k^2)*integral_0^t ds U(s)*integral_0^s dr U(r)
                  *[cos(2*k*(s-r))-cos(2*k*(t-r))].

The normalization is `F=2*k*|v|^2`. For example, a constant potential
with the same oscillator initial data gives exactly

    F(t) = 1-U/(k^2+U)*sin(sqrt(k^2+U)*t)^2,
    coefficient_U^2 F = sin(k*t)^2/k^4-t*sin(2*k*t)/(2*k^3).

This constant-potential identity is only a quadratic normalization
control; a nonzero constant potential is excluded from the endpoint
hypothesis used for the higher-order derivative-transfer proof.

Insert the Abel regulator `exp(-epsilon*k)`. For `0<a<=b`,

    integral_0^infinity exp(-epsilon*k)*(cos(2*k*a)-cos(2*k*b))/k dk
       = (1/2)*log((epsilon^2+4*b^2)/(epsilon^2+4*a^2)).

The regulated kernel is nonnegative and at most `log(b/a)`. Taking
`a=s-r`, `b=t-r` proves

    R_2[U](t) = integral_0^t ds U(s)*integral_0^s dr U(r)
                           *log((t-r)/(s-r)),
    R_2'[U](t) = integral_0^t ds U(s)*integral_0^s dr U(r)/(t-r).

The logarithmic kernel is locally integrable on the time triangle. Its
mass is `t^2/2`; the derivative kernel has mass `t`. The `s=t` boundary
of the original logarithmic kernel is zero, also before removing the
regulator. Both formulas therefore hold for continuous potentials.

The limiting operation is quantitative in `C1`. With `c=epsilon/2`, the
positive error in the first kernel is bounded by
`log(1+c^2/(s-r)^2)/2`. Its time-triangle integral is at most
`pi*epsilon*T/4`, using
`integral_0^infinity log(1+c^2/x^2) dx=pi*c`.
The derivative kernel error has triangle integral

    integral_0^t c^2/(b^2+c^2) db <= pi*epsilon/4.

After multiplication by `B^2`, these give the respective uniform errors.
Thus taking the time derivative and removing the regulator is justified,
not merely formal.

Polarizing the product of potentials gives

    ||R_2'[U]-R_2'[V]|| <= 2*B*T*||U-V|| <= 2*z*D.

## 4. All higher orders with one derivative only

Choose an arbitrary splitting frequency `K>0`. On `0<=k<=K`, use
`|sin(k*r)/k|<=r<=T`; on `k>=K`, use `|sin(k*r)/k|<=1/k`. The ordered
simplex volume and the product binomial sum give

    |w_n| <= (B*T^2)^n/n!                      (infrared),
    |w_n| <= (B*T/k)^n/n!                      (ultraviolet),
    |F_n| <= (2*B*T^2)^n/n!                    (infrared),
    |F_n| <= (2*B*T/k)^n/n!                    (ultraviolet).

For every `n>=3`, the full radial integral is absolutely convergent and
is bounded by `c_n*B^n`, where

    c_n_IR = (K^2/2)*(2*T^2)^n/n!,
    c_n_UV = K^2*(2*T/K)^n/[n!*(n-2)],
    c_n = c_n_IR+c_n_UV.

The ultraviolet integral for `n=2` diverges if one discards its
oscillatory cancellation. Section 3 is therefore essential, and the
numeric interface rejects `n<3` here.

Every term of `F_n'` contains one factor `U'` and `n-1` factors `U`.
Polarization or a two-step polynomial derivative gives

    ||R_n'[U]-R_n'[V]||
      <= c_n*[n*B^(n-1)*D+n*(n-1)*B^(n-2)*M*||U-V||]
      <= n^2*c_n*B^(n-1)*D.

This proof controls the actual mode integrals uniformly, not just a
first variation at `U=0`. The estimates permit all sums, derivatives,
and frequency integrals in `sum_{n>=3} R_n` to be interchanged by
absolute uniform convergence. Together with the separate Abel argument
for `R_2`, they define a continuous `C1` nonlinear functional on the
stated domain, with the physical agreement on the smooth past-flat
subset specified in section 1.

Take `K=1/T`, put `q=2*B*T^2=2*z`, and use, for `n>=3`,

    n^2/n! <= (3/2)/(n-3)!,
    3/2-n/[(n-1)*(n-2)] = (n-3)*(3*n-2)/[2*(n-1)*(n-2)] >= 0.

Then

    sum_{n>=3} n^2*c_n*B^(n-1) <= 18*z^2*exp(2*z).

The result is

    ||R'[U]-R'[V]|| <= [2*z+18*z^2*exp(2*z)]*D,
    ||R[U]-R[V]|| <= T*[2*z+18*z^2*exp(2*z)]*D.

The second estimate follows by integrating from the common zero initial
response. No restriction on the size of `z` is needed for these analytic
exponential bounds. The numeric calibration below explicitly restricts
`z<=1/4` to replace exponentials by simple rational upper bounds.

## 5. Shared nonzero prehistory and a short future interval

Suppose the two potentials coincide on `[0,t0]`, and write `T=t0+L`,
with `0<L<=T`. The same nonzero prepared history is retained; there is
no change to a new plane-wave state at `t0`. Set `W=U-V`. On the final
interval,

    |W(t)| <= (t-t0)*D,
    ||W||_1 <= L^2*D/2, ||W'||_1 <= L*D.

### Marked-factor simplex estimate

The factor carrying `W'` or `W` is supported in the future interval.
For a chain with n time variables, summing the possible positions of a
marked factor at time `s` gives the elementary volume identity

    sum_{j=1}^n s^(j-1)*(T-s)^(n-j)/[(j-1)!*(n-j)!]
      = T^(n-1)/(n-1)!.

It follows that a polynomial bound `c_n*B^n` becomes
`n*c_n*B^(n-1)*||Z||_1/T` when one potential factor is replaced by `Z`.
For two ordered chains in the product defining `F_n`, the same result
holds: the two insertion choices and the product binomial sum produce
`2^n/(n-1)!` in place of `2^n/n!`. A second differentiated potential
factor contributes at most its sup norm and the factor `n-1`.

Consequently, the derivative difference at order n is bounded by

    c_n*n*B^(n-1)*(L/T)*[1+(n-1)*L/(2*T)]*D
      <= (L/T)*n^2*c_n*B^(n-1)*D.

Thus the higher-order contribution is at most

    18*B^2*T^3*L*exp(2*B*T^2)*D.

### Quadratic shared-history improvement

Use the exact derivative kernel from section 3 and expand
`U(s)U(r)-V(s)V(r)=U(s)W(r)+W(s)V(r)`. At a future time
`t=t0+ell`, the two terms are bounded respectively by

    B*ell^2*D/2,
    B*D*integral_0^ell (ell-r)*log(t/r) dr
      = B*ell^2*[log(t/ell)/2+3/4]*D.

The total is at most

    B*L^2*[5/4+log(T/L)/2]*D.

Indeed the relevant expression is increasing as a function of `ell`
on `(0,T]`. This proof keeps the interaction of a future variation
with all of the old history; discarding that interaction would be
incorrect.

Combining the two pieces gives

    ||R'[U]-R'[V]|| <= C_shared*D,
    C_shared = B*L^2*[5/4+log(T/L)/2]
               +18*B^2*T^3*L*exp(2*B*T^2),
    ||R[U]-R[V]|| <= L*C_shared*D.

In particular `C_shared -> 0` as `L -> 0` with the complete original
history fixed. This is a genuine shortening estimate in the uniform
norm of the potential derivative.

## 6. Rational calibration and dimensions

For `0<=z<=1/4`, the geometric series gives

    exp(2*z) <= 1/(1-2*z) <= 2.

For `r>=1`, put `x=(r-1)/(r+1)`. With `N>=1`,

    log(r) <= 2*sum_{j=0}^{N-1} x^(2*j+1)/(2*j+1)
                  +2*x^(2*N+1)/[(2*N+1)*(1-x^2)].

This comes from the positive atanh series and bounding the remaining
denominators below by `2*N+1`. Every numeric bound returned by
`mode_lipschitz.py` is therefore rational; binary floats, booleans,
nonfinite values and unresolved symbolic inequalities are rejected.

The genuine calibration uses the actual dimensionless conformal
coordinate `eta_actual/eta_star`, not the background label `y`:

    delta=10^-14, T=3, L=1/4,
    M=2*delta*(61013499/8192).

The free starting slice for this calibration is `eta_actual=eta_star`
(`y=1`). The originally specified state propagates freely from `y=1/2`
to that slice, since the potential vanishes there, so only its harmless
overall phase origin changes. This is not a state reset at the plateau.
Here `T=3` is a uniform upper bound for the active history ending at any
observation in `2<=y<=3`. If its actual duration is `T_obs<3`, the
estimates apply on `[0,T_obs]`, using `T=3` only in upper bounds on
integral lengths and simplex volumes. They require no continuation of
the physical potential, or of the A.7 derivative estimate, beyond the
observed interval. The same interpretation permits `L` to be an upper
bound on the future interval length. Apply the exact exponential/log
estimates, which are increasing in these length bounds, before replacing
their elementary functions by the rational majorants.

A.7 derives the parenthesized coefficient for the prepared potential's
first actual-conformal derivative. The factor two puts that potential
inside, rather than on the boundary of, the derivative ball. With the
rational exponential bound and eight log-series terms, exact arithmetic
gives

    full-history C1 coefficient < 10^-8,
    shared-history C1 coefficient < 10^-10.

These coefficients multiply `D=||U'-V'||`. They are not SEE contraction
constants, stress-tensor bounds, or errors of order epsilon in a solution.

Under a change of time unit `t -> c*t`, the potential and its derivative
scale as `U -> U/c^2`, `M -> M/c^3`. Thus `z=M*T^3` and the reported
derivative Lipschitz coefficients are dimensionless. `R` scales as
inverse time squared, while `R'` and `D` scale as inverse time cubed.
Restoring the nonlinear contribution to `a^2 W` multiplies the estimates
by `hbar/(4*pi^2)`; comparing `W` itself on two metrics additionally
requires the elementary bounds on the two factors `a^-2`.

## 7. Relation to the integrated trace and remaining exclusions

The relevant source route is Meda--Pinamonti--Siemssen,
[2007.14665](https://arxiv.org/pdf/2007.14665), Proposition 3.2 and
Eq. (35), which separate the quadratic term and transfer a time
derivative to the potentials. Proposition 4.3 reduces an auxiliary
integrated-trace equation to an ordinary forced second-order equation.
Lemma 5.7 uses `X=a''/a` and a bounded inverse of the remaining linear
logarithmic operator. Our bound is a new explicit massless IR/UV
calculation of the nonlinear term, not an invocation of that paper's
full existence theorem or its unspecified constants.

In the massless minimal case `X=-U`, so the norm just bounded is exactly
the uniform norm of `X'` used there. No positive frequency gap is needed
for this estimate: the direct infrared Volterra bound is regular at
`k=0`, even if `U` changes sign. The paper's initial gapped reference
mode construction is a different decomposition, and its strict
positivity assumption must not be silently dropped when citing its
theorem. A positive initial plateau `U` can supply that gap, but the
present full-history estimate avoids introducing that reference state.

The elementary auxiliary-equation and metric bounds can also be
obtained without a frequency gap. For example, if two metrics have
identical `a,a'` at the future starting slice and solve
`a''+U*a=0`, then, writing `a_max` for a common bound and `B` for a
common potential bound,

    ||a_U-a_V|| <= (L^2/2)*a_max*exp(B*L^2/2)*||U-V||,
    ||a_U'-a_V'|| <= L*[a_max*||U-V||+B*||a_U-a_V||].

These follow by the scalar Volterra series and its derivative. They
retain possible nonzero initial curvature. The printed estimate for
`a'-a0'` in MPS Lemma A.2 appears to omit a term linear in `L` when
`X0!=0`; for instance `a''=X0*a`, `a0>0`, `a0'=0` gives
`a'-a0'=a0*X0*L+O(L^3)`. We use the direct inequalities above instead;
this observation is only a caution about that displayed estimate, not
a claim that the paper's existence result is false.

There are still four distinct requirements for an actual smooth SEE
application:

1. The local geometric/scheme terms and auxiliary trace solution must
   be bounded on a positive-metric ball, in addition to `R`.
2. The complete retarded map needs explicit self-mapping and
   contraction constants on the intended interval.
3. The initial energy constraint must hold exactly with the actual
   state and fixed classical-radiation amplitude. A small residual is
   not that constraint. Keeping the old metric fixed up to the exact
   evolution slice also fixes its initial RSET; it cannot correct a
   pre-existing constraint mismatch.
4. A `C1` potential fixed point gives only a `C3` metric without an
   additional smoothness/Hadamard compatibility argument. Changing
   metric jets while holding an old Cauchy covariance fixed is not a
   proof of that compatibility. A joint state/metric preparation would
   have to be controlled, not silently substituted.

No result here asserts that these four requirements are already met.
