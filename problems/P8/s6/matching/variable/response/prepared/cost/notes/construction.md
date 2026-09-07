# An exact loading source with an explicit Sobolev cost

This component concerns only the actual **punctured delta=0 limiting TT
operator** of S6.21, on `J=(a_*,b_*)=(-1/50,-1/100)`, with
`K=(tau kcom)^2` in `[1,4]`. It does not assign an action to the singular
delta=0 bounce. The physical source is the unchanged external g-metric
probe of S6.21/S6.23. All derivatives below are with respect to
`u=T/tau`, except the explicitly scaled coordinate `x=(u-a_*)/ell`,
`ell=1/100`.

## 1. Literal source operator and full endpoint normalization

Put `d=1+u²`. The coefficients in the unit-polarization TT action are

```
K_g=d^6/8, K_f=1/(2d^6), G_g=d²/8, G_f=1/(2d²),
U=4(1-u²)/[d^8(d^4-1)], P=1/U=d^8(d^4-1)/[4(1-u²)].
```

Both kinetic coefficients and U are positive on J. The physical equations
and probe convention are

```
(K_g g')'+G_g K g+U(g-f)=a³ sigma/4,
(K_f f')'+G_f K f+U(f-g)=0,
sigma=tau² Pi/M²,
S_probe=(1/2) integral dT a³ Pi gamma_g.
```

Consequently the first-order physical vector `X=(g,g',f,f')` has input
`B=(0,2,0,0)`: `a³/(4K_g)=2` exactly. The stress is a linear external
transverse-traceless spatial probe, not either internal scalar source.
Its zero time components, zero spatial trace and transverse polarization
give background covariant conservation for every C1 time profile.

Let `z` denote the **four-component Frobenius-normalized endpoint target**.
Its actual physical endpoint is

```
X(b_*)=C0(b_*) W_-(ell) Psi_-(ell;K) z.
```

Here C0 is the complete moving canonical-to-physical derivative map,
W is the regular-singular radial map, and Psi is the **full** wave
operator. None is replaced by its center value, and `Psi L v` is not
replaced by `L v`. The frozen S6.21 bound is `||Psi_-(ell)||<2`.

For the scaled physical endpoint `Z=diag(1,ell,1,ell)X`, direct entries of
`diag(1,ell,1,ell) C0 W` are bounded by

```
[1,0,1/5,0],
[ell,ell,101/1000,13/20],
[1,0,1/20,0],
[ell,ell,13/500,13/80].
```

Indeed the field map is `g=a_g l+b_g Q`, `f=a_g l+b_f Q`, where

```
a_g=2d³/sqrt(d^12+4),
b_g=-4/[d³ sqrt(d^12+4)], b_f=d^9/sqrt(d^12+4).
```

On `|u|<=1/50`, `|a_g|<1,|b_g|<2,|b_f|<1/2`; all three first
derivatives have absolute value below one. The latter follow from the
respective logarithmic-derivative bounds `6|u|,18|u|,18|u|`.
At the endpoint `sqrt(ell)=1/10`, and `mu=sqrt(39)/2<13/4`.
The sum of the squares of the displayed entry bounds is below four.
Thus

```
||Z(b_*)||_infinity <= ||Z(b_*)||_2 <4||z||_2.       (1)
```

This is where the true endpoint frame enters the quantitative cost.

## 2. Endpoint jets and exact degree-eleven construction

Let A_x be the scaled physical generator. It follows from the displayed
physical equations, not a frozen mass-only oscillator. In particular

```
A_x = [[0,1,0,0],
       [-ell²(K/d^4+U/K_g),-ell K_g'/K_g,ell² U/K_g,0],
       [0,0,0,1],
       [ell² U/K_f,0,-ell²(Kd^4+U/K_f),-ell K_f'/K_f]].
```

All its derivatives needed here are evaluated at `u=b_*`, which is
strictly punctured. Put `T0=I` and

```
T_(n+1)=sum_(j=0)^n binomial(n,j) (partial_x^j A_x)(b_*) T_(n-j).
```

Then `F_j=e_f^T T_j Z(b_*)`, `j=0,...,5`, are the exact target f jets
with respect to x. This recursion includes every background coefficient
derivative and the noncommuting matrix order. There is no substitution
of constant coefficients before differentiating.

Define f on the loading interval as the degree-eleven Hermite polynomial
whose jets0 through5 vanish at x=0 and equal F0 through F5 at x=1.
An explicit Bernstein representation has coefficients B0 through B5 zero
and

```
B_(11-j)=sum_(k=0)^j (-1)^k binomial(j,k) F_k/(11)_k,
j=0,...,5,
```

where `(11)_k` is the falling factorial. Set

```
g=f+P[(K_f f')'+G_f K f],
sigma=(4/a³)[(K_g g')'+G_g K g+U(g-f)].             (2)
```

The unsourced f equation is an identity, and the g equation has exactly
the physical normalization above. Extend the fields by zero before J
and by the full homogeneous target after J. Matching f jets0 through3
matches both physical fields and both velocities. Matching through5
additionally matches sigma and sigma' to zero at both endpoints.
The extended source is therefore C1 and belongs to `H0²(J)`: its
piecewise classical second derivative is square integrable and carries
no delta distribution. Direct source-boundary jet checks are included.

Its support may meet the **boundary** of J; it is not advertised as
C-infinity with support strictly inside J. This distinction matters for
attainment questions. Density plus a finite-dimensional endpoint
correction can compare the corresponding infima, but is not silently
used as an equality of source classes here.

## 3. Continuous polynomial bounds, not sampled controls

Writing `g=t0 f+t1 f'+t2 f''` in (2) gives an exact scalar operator

```
sigma=sum_(m=0)^4 c_m(u,K) f^(m)(u).
```

The spring cancels: all c_m have only powers of `1+u²` and `1-u²`
in their denominators. Their coefficient bounds may therefore use the
filled mathematical interval `|u|<=1/50`, without claiming the original
delta=0 action exists at zero.

For `j=0,1,2` the following integers bound
`ell^(-m) |partial_u^j c_m|`, uniformly for all u in J and `1<=K<=4`:

| j / m | 0 | 1 | 2 | 3 | 4 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 14 | 91 | 33896 | 21745 | 10427 |
| 1 | 3 | 4733 | 11539 | 1139718 | 1087204 |
| 2 | 153 | 1258 | 606335 | 1272779 | 56985890 |

Each entry follows by the exact numerator coefficient-l1 sum at
`|u|=1/50,|K|=4`, divided by a strictly positive factorwise denominator
lower bound, and then upward integer rounding. No decimal approximation
or coefficient sample is used.

For the Hermite polynomial, differentiating its Bernstein basis and
bounding the resulting coefficient rows gives

```
||partial_x^j f||_infinity
 <= B_j ||Z(b_*)||_infinity,
(B0,...,B6)=(2,17,171,2971,36005,499025,5020519).   (3)
```

The Bernstein basis is nonnegative and sums to one on `[0,1]`.
Every coefficient row is a polynomial in K; its absolute coefficient
sum on `|K|<=4` bounds the complete momentum interval. The independent
Fraction engine obtains the same rows by solving a six-by-six **power
basis** Hermite system and then converting differentiated power
polynomials to Bernstein form. It separately composes the literal
differential operators instead of copying the production c_m formulas.

Combining the tables gives the exact finite constants

```
||sigma||_infinity <=445826321 ||Z(b_*)||_infinity,
||sigma_uu||_infinity <=762987669721436 ||Z(b_*)||_infinity.
```

In the second line Leibniz is applied as
`sum_(m=0)^4 sum_(r=0)^2 binomial(2,r) C_(2-r,m) ell^(-r) B_(m+r)`.
Thus time derivatives and interval-volume factors are both retained.
Equation (1), `|J|=ell`, and elementary Lp inequalities prove

```
||sigma||_infinity <2e9 ||z||,
||sigma||_1 <2e7 ||z||,
||sigma||_2 <2e8 ||z||,
||sigma_uu||_2 <4e14 ||z||.                        (4)
```

For the zero target all quantities vanish; equivalently one may state
the four bounds non-strictly for every target. The strict inequalities
in (4) are for nonzero z. All constants are sufficient and unoptimized.

## 4. A full four-dimensional right inverse and a Gramian lower bound

Let C_K map the physical `L²(J,du)` source to the **Frobenius-normalized**
right endpoint, from zero initial data. Construction (2) defines a linear
map R_K into `H0²(J)` with

```
C_K R_K=I_4,  ||R_K||_(R4 -> L²) <=2e8.
```

For any w, duality gives
`||w||²=<C_K^*w,R_Kw> <=2e8 ||C_K^*w|| ||w||`.
Consequently the endpoint loading Gramian satisfies the entirely
analytic bound

```
G_K=C_K C_K^* >=1/(4e16) I_4,  1<=K<=4.           (5)
```

This is a lower bound on the smallest eigenvalue, not its optimal value.
A separately validated Gramian computation may sharpen it. The
right-inverse and derivative cost in (4) do not depend on such a
computation. A source optimal in L² need not itself lie in H0².

## 5. Lower costs for true regular-light targets

Let `z=L v` with L the inclusion of the two regular light coefficients.
The exact light columns of Psi coincide with the delta=0 prepared
analytic even/odd solutions of S6.23. Denote their physical g fields by
g_e and g_o, with origin light data `(1,0)` and `(0,1)` respectively.
The conserved canonical form uses `pi_l=l'-2omega Q,pi_Q=Q'`.
Its light block is the unit symplectic matrix and its cross block with
the two Frobenius heavy modes vanishes. To justify the latter, the heavy
Frobenius columns have `Q=O(r^(1/2)), Q'=O(r^(-1/2))`, with light
components tending to zero by their defining Volterra normalization.
The full light equation then gives `l'=O(r^(3/2)), l=O(r^(5/2))`:
both relative forcing terms are `O(r^(1/2))` and the homogeneous light
integration constants vanish. The regular columns have `Q=O(r³)` and
`Q'=O(r²)`, while `omega=O(r)`. Every term in the light-heavy canonical
pairing therefore tends to zero as r decreases to zero. Conservation
makes that cross pairing exactly zero on J. The weaker O(r) light bounds
from the initial Volterra estimate already suffice for this vanishing;
no symplectic orthogonality is assumed just from a basis label.
Therefore the light components
of the actual loading kernel are exactly

```
h1(u)=-a³ g_o(u)/2, h2(u)=a³ g_e(u)/2.            (6)
```

This uses the complete solutions and physical source projection, not
the bare common-field ansatz. Both light moment identities remain
necessary even when the two heavy endpoint coefficients are required
to vanish as well.

The S6.23 coefficient-l1 bounds have radius1/20, so they apply throughout
J, not only its smaller final physical slab. At delta0, put rho=2/5.
The minimum u degrees give

```
|l_e|<=1+(1/40)rho², |l_o|<=1/50+(1/800)rho³,
|q_e|<=(3/400)rho², |q_o|<=(1/100)rho,
D<1/300, a³<101/100.
```

Together with the physical-map bounds this yields
`|g_e|<101/100`, `|g_o|<21/1000`. Hence

```
|h1|<2121/200000, |h2|<10201/20000.
```

Cauchy-Schwarz now proves, for any exact source preparing L v,

```
||sigma||_2 >= max((2000000/2121)|v1|,
                  (200000/10201)|v2|),
||sigma||_2 >19||v||,  ||sigma||_2 >900|v1|         (7)
```

for the respective nonzero targets/components. These are necessary
lower bounds, not near-optimality statements. The stronger lower bound
may be used when a component vanishes; no positive bound is claimed for
the zero target.

## 6. Physical units and the remaining matching obligation

The norms in (4),(7) are in du. For example
`||Pi||_L²(dT)=M² tau^(-3/2)||sigma||_L²(du)` and
`||partial_T² Pi||_L²(dT)=M² tau^(-7/2)||sigma_uu||_L²(du)`.
Every source and response can be scaled
by an independent linear-probe amplitude. These identities do not bound
nonlinear backreaction or supply a cutoff.

For a regular target L v, using the same delta0 Hermite source for the
actual positive-delta family is covered by the pinned S6.23 estimate

```
endpoint error <=delta[8600||v||+126000000||sigma||_1]
              <=delta(8600+2520000000000000)||v||.
```

For example delta<=1e-20 gives error below `2.53e-5||v||` in the stated
fixed endpoint frame. This supplies a concrete preparation/source-norm
interface for one selected sector. It does not prove low temporal
frequency, arbitrary-source EFT locality, original DHOST-row matching,
vacuum/state selection, scalar/vector health, or the full adopted B
gate. The explicit second-derivative bound can instead be used in a
separately specified spectral-cost functional.
