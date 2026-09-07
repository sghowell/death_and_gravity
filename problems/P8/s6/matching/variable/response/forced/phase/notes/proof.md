# Fixed-pulse phase separation in the actual physical response

All variables, source normalization and four zero initial data are those
in [FORMULATION.md](../FORMULATION.md). This is a written analysis proof,
with exact replay of its algebra and rational margins, not a formalized
proof or a numerical delta scan. Fix r, K and the source profile before
taking delta to zero. We write g_delta for gamma_g,delta(r).

## 1. A comparison response with all omitted terms bounded

The frozen S6.26 theorem gives the exact light inverse G_L, the universal
retarded heavy inverse G0 and

```
l0=G_L(jL sigma),
gapp=ag*l0+bg*G0(jH sigma-Bop*l0),
||g-gapp||infinity<=1000*r^4.
```

Its proof also bounds `||G0||<=101*r^2/75`, `|bg|<2`,
`N(l0)<5*r^2`, and `||Bop*l0||<5*r^2*(10+60*r^2)`, on the
whole interval and whole K band. Define the scalar endpoint comparison

```
H_delta=ag_delta(r)*l0_delta(r)
       +bg_delta(r)*integral_{-r}^{-r/2}
          G0_delta(r,s)*jH_delta(s)*sigma(s) ds.
```

Then, including both the full feedback and the bounded pole remainder,

```
|g_delta-H_delta|
 <=[1000+(202/15)*(10+60*r^2)]*r^4
 <1200*r^4.                                           (1)
```

Nothing in this estimate drops the source or replaces the actual physical
map. The ordinary light coefficients A, jL and ag are smooth in delta
down to zero on this fixed compact interval. Their Volterra equations
therefore imply `ag_delta(r)*l0_delta(r)->L_K`, a finite real number.
For example, subtract the two light Volterra equations; uniform coefficient
convergence and their uniformly convergent factorial majorants give this
conclusion for every bounded fixed sigma. No heavy singular coefficient
occurs in that light equation.

## 2. The finite-clock cross-bounce asymptotic

Put `mu=sqrt(39)/2`, `rho_delta(u)=sqrt(u^2+delta/8)` and
`z_delta(u)=asinh(sqrt(8)*u/sqrt(delta))`. The S6.26 kernel is

```
G0_delta(r,s)=sqrt(rho_delta(r)*rho_delta(s))/mu
             *Im[fplus(z_delta(r))*conj(fplus(z_delta(s)))].
```

The single Gauss solution is
`fplus(z)=exp(i*mu*z)*2F1(-1/2,3/2;1-i*mu;(1-tanh(z))/2)`.
The standard connection between arguments zero and one gives

```
fplus(z)=exp(i*mu*z)*(1+o(1)),                         z->+infinity,
fplus(z)=A*exp(i*mu*z)*(1+o(1))
        +B*exp(-i*mu*z)*(1+o(1)),                     z->-infinity,

A=Gamma(1-i*mu)*Gamma(-i*mu)
  /[Gamma(3/2-i*mu)*Gamma(-1/2-i*mu)],
B=Gamma(1-i*mu)*Gamma(i*mu)/[Gamma(-1/2)*Gamma(3/2)].
```

Here c-a-b=-i*mu is nonintegral; the connection uses no pair of solutions
at the hypergeometric singularity infinity, where a-b is integral. Both
local Gauss series at zero have nonexceptional denominator parameters.
To see this for both left terms, put x=(1-tanh(z))/2 and y=1-x.
Euler's identity gives the exact connection

```
fplus(z)=A*exp(i*mu*z)*F(-1/2,3/2;1+i*mu;y)
        +B*exp(-i*mu*z)*F(-1/2,3/2;1-i*mu;y).
```

In the second term the original Gauss factor acquires x^(i*mu);
the identity y/x=exp(2z) combines this with y^(-i*mu) into the displayed
plane wave. The coefficient-ratio bound for either resulting series is

```
|(n-1/2)*(n+3/2)/[(n+1)*(n+1+-i*mu)]|<=1, n>=0.
```

It follows by expanding the squared-denominator minus squared-numerator
as a polynomial with positive coefficients, checked independently. Thus
`|F(q)-1|<=q/(1-q)` and `|d_z F(q)|<=2*q/(1-q)^2` when the small
argument is x or y. It is the physical-z derivative remainder that tends
to zero; the Gauss-argument derivative F_q is generally only bounded.
At the right endpoint and left source interval all small arguments obey
`q<=delta/(8*r^2)<=1/800`. These function and physical-z derivative
remainders therefore vanish uniformly. The real ODE potential
makes the Wronskian constant. Evaluating it on either side gives
`|A|^2-|B|^2=1`; in particular `|A|>=1`. No fitted Gamma modulus is used.
These are auxiliary ODE basis statements, not quantum state assumptions.

For v in [r/2,r], elementary real-log asymptotics give, uniformly in v,

```
z_delta(r)-z_delta(-v)=log(32*r*v/delta)+o(1),
z_delta(r)+z_delta(-v)=log(r/v)+o(1),
sqrt(rho_delta(r)*rho_delta(-v))=sqrt(r*v)+o(1).
```

All inputs are bounded away from zero at fixed r. More explicitly,
`0<=asinh(w)-log(2*w)<=1/(4*w^2)` for w>0 follows from the square-root
and logarithm elementary bounds. Setting w=sqrt(8)*v/sqrt(delta) gives
a clock remainder at most delta/(32*v^2). The positive fourth-root
factors in the amplitude also converge uniformly, by concavity. These
estimates and the preceding Gauss bounds give a uniform error tending
to zero, without exchanging a fixed-r limit with r=0. Consequently

```
G0_delta(r,-v)=sqrt(r*v)/mu *Im[
   conj(A)*exp(i*mu*log(32*r*v/delta))
  +conj(B)*exp(i*mu*log(r/v))] +o(1),                 (2)
```

with a uniform remainder. This uniformity justifies integrating (2)
against the fixed bounded pulse. The actual weights jH_delta converge
uniformly on the source interval, and bg_delta(r) converges too.

## 3. A nonzero oscillation for every admitted pulse

Use the actual limiting weights

```
j0(v)=-2*(1+v^2)^3/sqrt((1+v^2)^12+4),
bg0(r)=-4/[(1+r^2)^3*sqrt((1+r^2)^12+4)],
I_sigma=integral_{r/2}^r sqrt(v)*j0(v)*sigma(-v)
                       *exp(i*mu*log(v/r)) dv.
```

Equations (1)-(2) give

```
H_delta = D_K + Im[C_sigma*exp(i*mu*log(32*r^2/delta))]+o(1),
C_sigma=bg0(r)*sqrt(r)*conj(A)*I_sigma/mu,             (3)
```

where D_K is real and delta-independent. In particular, the B term
contributes only to D_K. Small reflection is not needed for this proof;
the transmission amplitude is nonzero even if B were zero.

First replace j0 by -2/sqrt(5) and sigma by 1 on [r/2,r]. The exact
dimensionless moment is

```
integral_{1/2}^1 sqrt(x)*exp(i*mu*log(x)) dx
 =[1-2^(-3/2-i*mu)]/(3/2+i*mu).
```

Its modulus is at least `(1-2^(-3/2))/sqrt(12)>5/28`, because
`2^(-3/2)<3/8` and `sqrt(12)<7/2`. Pulse smoothing changes this normalized
moment by at most 1/100: use sqrt(v)<=sqrt(r) and the declared L1 deficit.

The source-weight correction is also uniformly small. With t=v^2 and
q=(1+t)^6,

```
(j0(v)/(-2/sqrt(5)))^2=5*q/(q^2+4),
5*q-q^2-4=(q-1)*(4-q)>=0.
```

On our domain q<4. The positive source ratio is therefore at least 1;
it is at most `(1+t)^3<=1+4*r^2`, by bounding its denominator below.
Since sigma<=1 and the interval length is less than r, its additional
moment error is at most `4*r^2` in the same units. Hence

```
|I_sigma| > (2/sqrt(5))*r^(3/2)*(5/28-1/100-4*r^2).
```

Likewise `(1+r^2)^12+4<=5*(1+r^2)^12` gives
`|bg0(r)|>=(4/sqrt(5))*(1+r^2)^-9
 >=(4/sqrt(5))*(1-9*r^2)`.
The final inequality follows either by differentiation or by multiplying
through and checking the nonnegative binomial coefficients, replayed
independently. With `mu<13/4` and `|A|>=1`,

```
|C_sigma|/r^2
 >(32/65)*(1-9*r^2)*(5/28-1/100-4*r^2)
 >2/25.                                             (4)
```

Both positive factors in the middle decrease with r^2. Their exact
values at r=1/1000 prove the strict lower bound for every smaller r;
this is a continuous inequality, not a sample-grid inference.

## 4. Separation of the actual coupled responses

C_sigma is nonzero. Choose its argument theta. The positive sequences

```
delta_n^+=32*r^2*exp[-(pi/2-theta+2*pi*n)/mu],
delta_n^-=32*r^2*exp[-(3*pi/2-theta+2*pi*n)/mu]
```

tend to zero, and eventually satisfy the original delta domain. Equation
(3) approaches D_K+|C_sigma| and D_K-|C_sigma| on these sequences.
We do not claim the actual responses themselves converge on either
sequence: their residuals in (1) need only be bounded. Those bounds suffice:

```
limsup g_delta-liminf g_delta
 >=2*|C_sigma|-2400*r^4
 >(4/25-2400*r^2)*r^2
 >3*r^2/20.                                         (5)
```

Boundedness follows from (1)-(3), or directly from the S6.26 operator
norms. Thus the limsup and liminf are finite. The source, measurement
time, physical map prescription and K are fixed in this comparison;
only the parameter delta of the recorded family changes.

## 5. Smooth sources and the exact exclusion boundary

Let eta(x)=exp(-1/x) for x>0 and zero otherwise, and
S(x)=eta(x)/(eta(x)+eta(1-x)). Set w=r/800 and

```
sigma(-v)=S((v-r/2-w)/w)*S((r-w-v)/w).
```

This is smooth with compact support strictly inside (r/2,r) in the v
coordinate, lies in [0,1], and equals 1 except within strips of total
length at most 4w=r/200. It satisfies the required deficit and is
independent of delta. Its extension by zero gives a smooth conserved
external TT probe with a flat-past preparation. This establishes smooth
examples without interpreting a distributional impulse as a low-band
source. Arbitrarily small positive source rescalings preserve a
correspondingly rescaled separation.

The excluded statement is convergence, for this class of fixed pulses,
of the actual endpoint response to a delta-independent limiting response.
It does not exclude an explicitly phase-dependent effective model or
general source classes with different histories or initial correlations.
It does not establish a rolling gap, a Wilsonian cutoff or a
parametrically small temporal-source hierarchy. The full-parent center
algebraic mass alone therefore does not justify this particular
fixed-source limit. No original P8, C/D or UV closure follows.
