# P8(a) A.8 — Actual reference-stress and SEE-residual bounds

Adopted 2026-09-06. The immutable immediate dependency is A.7, certificate
SHA256 `cc4ec02bcb23e0ef01d0fb58ad4b95a3bb602a5bbbfc6871d044b2a09a78e963`.
No earlier checkpoint is changed.

## Result and fixed physical model

On the **same** prepared metric and state as A.7, this checkpoint bounds
the actual renormalized density, pressure and EED against the positive A.3
radiation reference at the same proper time. It then gives a quantitative
finite-amplitude residual for the actual Einstein-plus-radiation equation,
including physical source amplitude epsilon=1 in an explicit regime.
Neither the metric nor a nearby metric is asserted to solve that equation.

Keep one real free massless minimally coupled scalar and the positive,
zero-mean homogeneous isotropic quasifree Hadamard state transported from
the unchanged radiation past. The past Cauchy slice may be y=1/2. All
potential jets vanish before preparation at y=1. Do not replace this state
by a Born covariance or by instantaneous flat modes on a curved metric.

Use the A-track FK convention `+---`, `R=+6*(Hdot+2H²)`. The source GS
Box equals FK Box; `g_FK=-g_GS`, `R_ab,FK=-R_ab,GS`, `R_FK=R_GS`.
Physical density and pressure agree, whereas the trace and finite-tensor
dictionary is `Theta_FK=-Trace_GS`, `I_GS=-I_FK`. The physical finite
coefficient is `gamma=-(c3_GS+c4_GS/3)` on FLRW.

The numerical theorem uses the already named raw-Hadamard prescription

    lambda=2*sqrt(2)*A*eta_star²=4*sqrt(2)*t_star, gamma=0,
    t_star=A*eta_star²/2.

There is no additional finite Wick-square `alpha*R` shift. The generic
symbolic reconstruction retains gamma and the homogeneous integration
constant explicitly, but the numerical calibration does not cover another
choice without adding its local contribution. Absolute stress is not
scheme-independent.

## Unchanged metric, amplitude and exact clock

Let A, eta_star, kappa and hbar be positive. Keep
`d=kappa*hbar/(46080*pi²)` and `epsilon>=0`, with

    y=sqrt(2*t/A)/eta_star,
    delta=16*epsilon*d/(A²*eta_star⁴)=4*epsilon*d/t_star²,
    chi_tilde(y)=sigma(y-1)*sigma(4-y), p(y)=chi_tilde(y)/y⁴,
    f(y)=1-delta*p(y), a=A*eta_star*y*f(y)^(1/4),

where sigma is exactly A.6/A.7's C-infinity switch. In proper time this is
the original `a_epsilon=sqrt(2*A*t)*(1-4*epsilon*d*chi(t)/t²)^(1/4)`.

The actual conformal derivative is `d_eta=f^(1/4)*d_y/eta_star`, and its
time differences are exact integrals of `eta_star*f^(-1/4)`. No inverse
clock or epsilon-linearized metric is used. Keep the derived A.7 domain

    0<=delta<=delta_bar=2048/190897521.

The observation envelope is `2<=y<=3`; compact interior targets retain
the original near-target plateau condition. All history bounds start at
y=1, not at the earlier Cauchy slice, because the potential and remainder
data are identically zero before the preparation begins.

## Complete trace/conservation reconstruction

Write `h=a'/a`, `U=-a''/a`, with primes now in actual conformal time. A.7
gives the exact Wick square

    W=<phi²>_ren=hbar/(4*pi²*a²)*(Rmode-K/2),
    K=integral U'(s)*log((eta-s)/T) ds
          +U(eta)*(log(sqrt(2)*a*T/lambda)+5/6).

The raw conserved trace is

    Theta=-Box W/2-hbar*v1/(4*pi²)-6*gamma*Box R,
    v1=R²/288+(Riem²-Ric²)/720-Box R/120.

The Box-R term is absorbed algebraically, without changing prescription,
by using `W_eff=W-hbar*R/(240*pi²)`. Since `a²*R=-6*U`, define

    S=Rmode-K/2+U/10,
    J=integral U'*S d_eta, L=integral h*U² d_eta,

with both histories zero in the actual common past. In units
`hbar/(pi²*a⁴)` the gamma=0 components are exactly

    rho=h⁴/960+[-h*S'+(h²-U)*S+J]/8-L/32,
    p=(5*h⁴+4*h²*U)/2880
          +[S''-3*h*S'+(3*h²+U)*S+J]/24+(U²-L)/96,
    E=(3*h⁴+2*h²*U)/960
          +[S''-4*h*S'+4*h²*S+2*J]/16+(U²-2*L)/64.

For generic gamma add `gamma*I_R²` to the physical components. A generic
conserved reconstruction also permits the physical homogeneous term
`(C/a⁴,C/(3*a⁴),C/a⁴)`. The specified actual past density fixes C=0:
there S and all its jets, J, L and I vanish, and the remaining Euler term
is precisely the pinned A.3 density. Conservation alone does not fix C.

## Quantitative actual-stress comparison

A.7 already derives every potential norm needed here. The proof shows
`|K^(j)|<=delta*L_j/eta_star^(j+2)` for j=0,1,2, with

    L0=3*b1+3*b0,
    L1=3*b2+3*b1+2*b0,
    L2=3*b3+3*b2+4*b1+(4+delta_bar*b0)*b0.

Combining these with the A.7 nonlinear remainder and the exact anomaly
gives, componentwise throughout the envelope,

    |T_actual-T_ref| <= hbar*(delta*C1+delta²*C2)/(pi²*A⁴*eta_star⁸).

The following convenient rounded coefficients strictly exceed the exact
rational coefficients pinned in the report:

| Component | C1 | C2 |
| --- | ---: | ---: |
| Density | 14000 | 800000000000000 |
| Pressure | 236000 | 900000000000000 |
| EED | 360000 | 1700000000000000 |

The reference is A.3 at the same proper time, namely
`hbar/(pi²*A⁴*eta_star⁸*y⁸)` times `(1/960,1/576,1/320)`.
For every `0<=delta<=10^-14`, exact rational comparisons give less than
one-percent error relative to each positive reference component. Thus the
actual **specified transported state's** density, pressure and EED are
positive on this target. There is no positivity assertion for arbitrary
Hadamard states.

## Quantitative finite-amplitude actual SEE residual

Use the existing source model: zero cosmological term, separately conserved
ordinary radiation `rho_rad=3*A²/(kappa*a⁴)`, `p_rad=rho_rad/3`, and
the actual quantum source multiplied by epsilon. No additional explicit
curvature-squared gravitational sources are included in this calibration.
If they are retained, their contributions require additional bounds and
cannot be silently absorbed into the reported constants.

Define the physical-sign residuals

    D_rho=3*H²-kappa*(rho_rad+epsilon*rho_actual),
    D_p=-2*Hdot-3*H²-kappa*(p_rad+epsilon*p_actual),
    D_E=(D_rho+3*D_p)/2.

The corresponding FK equation components have opposite signs and the same
absolute bounds. On the target plateau, put z=delta/y⁴. A.5's exact frozen
defects and the actual stress bounds imply

    |D_i| <= delta²/(A²*eta_star⁴)
         *[c_i/1922+2880*(C1_i+delta*C2_i)], c_i=(3,5,9).

The denominator is literally `1922=2*31²`. The proof derives it from
`z<=1/32`, `y¹²>=4096`, and `(2-z)/(1-z)²<=2048/961`.
The coupling restoration is exact:
`epsilon*kappa*hbar/pi²=2880*delta*A²*eta_star⁴`.

At `delta<=10^-14` the residuals are each below `10^-17` of the respective
positive zeroth-order classical Einstein/radiation component on the envelope.
These comparison scales are known reference quantities, not unknown actual
residuals or small quantum quantities used as denominators. Physical
epsilon=1 is included whenever `t_star>=2*10^7*sqrt(d)`.

## Scope boundary

This closes the numerical actual-response input for the stated fixed
prepared testbed and source model. It does not solve SEE, prove closeness
to a solution or stability, impose self-consistent initial data, or control
the residual over the off-shell preparation interval. It does not change
the metric to A.6's second-order corrected family. No full two-point
function QSEI estimate, arbitrary-Hadamard positivity, all-geodesic focusing,
new singularity theorem or full P8(a) completion is claimed.
