# Full spatial flat covariant/time-subtraction conversion

## Same parent, explicit benchmark

Keep the canonical ordinary Proca sector with m=1000, kappa=10^800
and the S6.199 full flat stress cut. The actual CD prepared state and
original finite prescription are unchanged. The benchmark uses arbitrary
compact smooth symmetric spatial metric tests; its nonlocal cut does
not fix the metric-parametrization-dependent local contacts.

Let p be the spatial momentum, q=|p|^2, w=(omega+i0)^2 and z=w-q.
For spatial indices define A(z)=z I+p p^T,

    Q0_ijkl=A_ij A_kl/3,
    Q2_ijkl=(A_ik A_jl+A_il A_jk)/2-Q0_ijkl,
    N_i(z,p)=z Q_i(z,p).

Q_i=z^2 Pi_i for the restricted full covariant spin projectors.
The polynomial numerator N has no light-cone pole.

## Exact conversion

Define the explicit covariant nonlocal representative

    F(z,p)=sum_i N_i(z,p) integral rho_i(s)/(s^3(s-z)) ds,

with both S6.199 densities and s from 4m^2 to infinity. This chooses
a nonlocal representative, not the physical finite local polynomial.
Write N_i(w-q,p)=sum_(k=0)^3 n_ik(p)w^k and

    J_in(q)=integral rho_i(s)/(s^3(s+q)^(n+1)) ds.
    A_r(p)=sum_i sum_(k=0)^r n_ik(p)J_i,r-k(q), r=0,1,2.

Then exactly

    F(w-q,p)=A0+w A1+w^2 A2+B6(w,p),
    B6=w^3 integral sum_i rho_i(s)Q_i(s,p)/s^2
                  /((s+q)^3(s+q-w)) ds.

B6 is the full flat retarded sixth-time-derivative subtracted sine
bulk. The conversion acts as A0 Gamma-A1 Gamma''+A2 Gamma''''.
Each A_r is finite and time local but is generally not a finite
spatial differential operator. Both spin sectors and the entire
spatial tensor, not just a TT projection, enter the identity.

## Quantitative continuum result

At every spatial momentum,

    ||A_r|| <= (q+m^2)^3/(80pi^2 m^(2r+2)),
    ||A_r-A_r,Lambda|| <= (q+m^2)^3/(9pi^2 m^(2r) Lambda),

for the computational squared-spectral-mass limit Lambda>=4m^2.

On the unit slab let

    Z(Gamma)^2=sum_(r=0)^2||(1-Delta)^3 partial_t^(2r)Gamma||L2^2,
    U(Gamma)^2=Z(Gamma)^2+||partial_t^6 Gamma||L2^2.

The conversion has weak bound 1e10||D||L2 Z and tail
1e17||D||L2 Z/Lambda. The absolute sine-bulk bound is
1e-5||D||L2||Gamma^(6)||L2, with tail
||D||L2||Gamma^(6)||L2/(100sqrt(Lambda)).
Consequently the complete specified F has bound

    2e10 ||D||L2 U(Gamma),

and full spectral-tail error

    1e14 ||D||L2 U(Gamma)/sqrt(Lambda).

Both external-metric canonical factors give 8e-790 and
4e-786/sqrt(Lambda) in the same norms.

## Unclosed boundary

These are derivative-losing weak flat nonlocal-response estimates,
not a full constraint-reduced mixed norm or same-space contraction.
No physical cutoff, finite local polynomial, curved odd-endpoint
discard, finite-regulator identification, actual CD spatial matching,
response inverse, nonlinear interacting background, stability or
original V/G/B closure follows. Scoped P8(a) and A.20-A.23 are unchanged.
