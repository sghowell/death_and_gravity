# S5.3a — Exact mode frequencies and uniform local indicators

Date: 2026-09-04. This is an independently checked local algebraic and
inequality calculation, not external peer review, a global WKB theorem or
an interacting strong-coupling certificate.

## 1. Include the expanding measure before reading a frequency

For a single diagonal scalar/tensor mode, write the quadratic Hamiltonian
density as `A p^2/2+B p Q+C Q^2/2`. Here the canonical momentum density is
`a^3 p`, and the action is `integral dt a^3[p Qdot-h]`. Eliminating p gives

```
L=a^3[K(Qdot-BQ)^2-CQ^2/2],   K=1/(2A).
```

The exact oscillator coordinate and its logarithmic derivative are

```
Y=a^(3/2) sqrt(2K) Q,
f=3H/2+(dot K)/(2K),   b=f+B.
```

After substitution and one integration by parts,

```
L=Ydot^2/2-Omega^2 Y^2/2 + time boundary,
Omega^2=A C-b^2-dot b.
```

The `3H/2` contribution must not be dropped. A generic symbolic quadratic
completion and independent boundary identity verify this formula. The
scalar calculation is also checked against direct differentiation in the
original uncompactified cosmic time at fixed comoving momentum, without
using the compact derivative operator.

## 2. Apply the exact S5.2 Hamiltonians

In the earlier fixed local units, `H=Theta=2x`, `J=2P` and
`Lambda=1-2(1-x^2)^3`. The scalar coefficients are

```
unitary:
  K=J/Theta^2,
  B=-Theta Lambda q/J,
  C=2(Lambda^2 q^2/J-q);
gamma:
  K=qJ/(q Lambda^2-J),
  B=Theta Lambda q/J-H,
  C=2 Theta^2 q^2/J.
```

The gamma `-H` term is the time-dependent canonical swap contribution.
For a kinetically normalized tensor polarization choose `K=1/2`, `B=0`,
`C=q`; its frequency is independent of that constant normalization choice.

Let `Omega^2 ell^2=q+mu(x,q)`. To differentiate compact coefficients use

```
D_w=(1-x^2)partial_x-2xq partial_q-2wx,
D_w f=ell^(2w+1) partial_t[f/ell^(2w)] at the evaluation point.
```

The scalar `dot K/(2K)` in `f` uses weight 0 in unitary and weight 1 in
gamma, as proved in the physical-reduction note. The dimension-one `b`
uses weight 1/2. This holds `ell` fixed in the patch rather than treating
the compact parametrization as a new global time coordinate.

For unitary the q terms cancel to exactly the old luminal principal q;
the remaining mass is `mu=-f^2-D_(1/2)f`. The exact tensor expression is
particularly simple:

```
mu_tensor=-3-3x^2.
```

The report stores the full scalar rational functions, not a high-q
truncation. At the literal bounce the gamma oscillator is finite:

```
Omega_gamma^2 ell^2=(q^3-16q^2+26q+4)/(q-2)^2
                 =q-12+O(1/q).
```

At either compact tail endpoint the unitary and tensor frequencies tend
to `q-6`. The unitary oscillator is singular at x=0 and is not used there.
The gamma velocity chart's denominator q=2 at the bounce is a chart issue,
not an inferred strong-coupling cutoff. None of these mass corrections
reverses the already certified *principal* luminality statement.

## 3. Compactify high frequency and bound the mass derivatives

Use `r=1/q`, `0<=r<=1/68`. The derivative becomes
`D_w=(1-x^2)partial_x+2xr partial_r-2wx`.
For each chart compute exactly

```
mu,    mu_1=D_1 mu,    mu_2=D_(3/2) mu_1.
```

These are the compact versions of the physical mass correction and its
first two cosmic-time derivatives. They have finite r=0 limits. The
following denominator and numerator bounds cover every real time:

- Exterior unitary domain `x^2>=1/65`: `|x|>1/9`, `P>1/4`, `|x|<=1`.
- Core gamma domain `x^2<=1/17`: `|x|<1/4`, `P>1/4`, `P<=1`,
  `Lambda^2>1/4`; hence `R=Lambda^2-2Pr>15/68` for `r<=1/68`.
- Tensor domain: all `|x|<=1`; its corrections are polynomials.

Factor each rational denominator exactly. It must be a nonzero constant
times powers of the proved factors `x,P` or `P,R` (or a constant for
tensors); any new factor raises an error. Majorize the numerator by the
sum of the absolute values of its power coefficients times the domain's
`|x|` and r upper bounds. Divide by the product of the denominator lower
bounds. This is a direct exact-rational triangle inequality over the entire
domain, not a grid sample. The previously certified P and Lambda bounds
are reused with their report/source pins.

The report records all rational majorants `C0,C1,C2` with
`|mu|<=C0`, `|mu_1|<=C1`, `|mu_2|<=C2`. The largest values come from
the deliberately crude unitary denominator weakening:

```
C0_unitary=1547424,
C1_unitary=7685736192,
C2_unitary=53229160793088.
```

These are sufficient bounds, not extrema or measured physical scales.
The gamma majorants are below 8619, 1319608 and 412614378 respectively;
the exact rational values, not these rounded-up summaries, are stored.

## 4. Explicit first/second adiabatic indicators

Let `V=Omega^2 ell^2=q+mu` and suppose
`q>=max(2C0,C1,C2,10^8)` in the selected chart. Then

```
q/2<=V<=3q/2,
|ell^3 partial_t Omega^2|=|-4xq+mu_1|<=5q,
|ell^4 partial_t^2 Omega^2|=|(-4+24x^2)q+mu_2|<=21q.
```

The exact chain rule for the positive root `Omega=sqrt(Omega^2)` yields

```
|dot Omega|/Omega^2 <= 5 sqrt(2)/sqrt(q),
|ddot Omega|/Omega^3 <= 92/q.
```

In the second bound the two terms are at most `42/q` and `50/q`.
The report selects the first power of ten above all required majorants,
giving **q>=10^14** uniformly across the covering charts and tensors.
Consequently the first indicator is below `7.1e-7` and the second at most
`9.2e-13`. These ratios and the positivity of V are proved inequalities;
the enormous sufficient threshold is not claimed to be optimal.

The scale comparison is local: `k_physical=sqrt(q)/ell`, while the prior
curvature scale is at most `sqrt(6)/ell`. Thus this diagnostic band is well
above curvature. It is **not yet shown to be below an interaction cutoff**.
A fixed comoving mode redshifts and need not stay in the band for all time.
The statement is that this local-frequency diagnostic band is available
at each background time with uniform constants, including both tails.

## 5. Remaining work is substantive

Small pointwise indicators alone do not specify a global vacuum, a uniform
solution error over an arbitrarily long interval, or a scattering matrix
on this time-dependent spacetime. Interaction coefficients and the temporal
window of a local amplitude still need their own variation/error treatment.
The S5.2 vertices have not yet been combined into physical exchange/contact
amplitudes; a nonempty adiabatic **and weakly coupled** window has not been
proved. The gamma/zero-spatial-momentum and massless-pole issues must be
handled with explicit frequency/angular/IR assumptions in that calculation.
No result here addresses loops, all orders, M1, nonlinear stability or UV
completion.
