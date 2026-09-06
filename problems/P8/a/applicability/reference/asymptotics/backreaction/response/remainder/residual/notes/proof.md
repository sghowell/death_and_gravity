# Exact reconstruction and finite-residual proof

All metric, history and quantum quantities refer to the actual A.7 clock
and state. A.7 supplies the full all-frequency C2 Wick remainder with its
derived potential-jet bounds; it is replayed rather than modified.

## 1. Anomaly decomposition and the homogeneous solution

In FK conventions put h=a'/a, U=-a''/a and H=h/a. The exact FLRW identity
`Riem²-Ric²=-12*H²*(Hdot+H²)` gives, in the named gamma=0 prescription,

    Theta=-Box(W-hbar*R/(240*pi²))/2
                -hbar*R²/(1152*pi²)
                +hbar*H²*(Hdot+H²)/(240*pi²).

Since `a²R=-6U`, the first term is
`-hbar*Box(S/a²)/(8*pi²)` with `S=Rmode-K/2+U/10`. This is an exact
algebraic rearrangement of the raw trace, not a change of finite Wick
renormalization. Omitting U/10 changes the trace on generic curved metrics.

The last trace term has the conserved density
`rho_Euler=hbar*H⁴/(960*pi²)`, because
`(d_t+4H)H⁴=4H³*(Hdot+H²)`. The middle trace term has

    rho_R²=-hbar/(32*pi²*a⁴)*L, L=integral h*U² d_eta,

because `R²=36U²/a⁴` and the conservation equation is
`(a⁴rho)'=a⁴h*Theta`. A.7's exact trace/conservation identity, now applied
linearly to S instead of Rmode, gives

    rho_S=hbar/(8*pi²*a⁴)*[-h*S'+(h²-U)*S+J], J=integral U'*S d_eta.

Combining these densities, then using `p=(rho-Theta)/3` and
`E=rho-Theta/2`, gives the complete displayed reconstruction. Direct
covariant symbolic tests and independent Laurent-polynomial controls check
the trace and conservation of the combined expression. Neither history may
be dropped: their derivatives cancel distinct local conservation terms.

The generic homogeneous density C/a⁴ is conserved and traceless, so trace
and conservation do not select it. In the actual common past, U, S and
all their required derivatives, J, L and I_R² vanish. The reconstructed
Euler density equals A.3's actual reference density there. Thus the
specified common past data, not a free choice made after integration,
fix C=0. A test deliberately changes C: conservation survives but the
physical initial density no longer agrees.

The generic symbolic interface also adds the finite physical gamma*I_R².
Its trace is `-6*gamma*Box R` and it is conserved. The numerical result
sets gamma=0 as part of the named prescription; the identity does not
make an arbitrary finite coefficient disappear.

## 2. Exact-clock logarithmic norm estimate

Choose the auxiliary span T=3*eta_star. The active history begins at y=1:
all potential/remainder data vanish between the earlier Cauchy slice and
this boundary. A.7 proves the active actual-conformal duration is <=T,
`|h|<=2/eta_star` throughout that history, and

    |U^(j)|<=delta*b_j/eta_star^(j+2),
    |h'|=|-U-h²|<=(4+delta_bar*b0)/eta_star².

The h cap is valid on the full active history, not only at observations:
`eta_star*h=f^(1/4)/y-delta*p'*f^(-3/4)/4`, y>=1, together with the
p-derivative amplitude cap proves it.

For the named physical length the local factor in K is

    ell_local=log((3/2)*y*f^(1/4))+5/6.

On 1<=y<=3, 1/2<=f<=1. Its log argument is greater than one because
`(3/2)^4/2>1`, and at most 9/2. The positive exponential series gives
`exp(2)>1+2+2=5>9/2`, so `0<log(argument)<2` and
`|ell_local|<17/6<3`. All displayed margins are rational. The length
is constant in time; `ell_local'=h`, `ell_local''=h'`.

At duration D<=T,

    integral_0^D |log(x/T)| dx = D*(1+log(T/D)) <= T.

The last expression increases with D on [0,T], has value T at D=T and
continuous value zero at D=0. Write the retarded integral with x=eta-s.
Its first two derivatives differentiate U', not the logarithmic singularity.
The moving-upper-limit terms vanish because U' and U'' at the past
boundary are zero. Smoothness and the integrable logarithm justify these
differentiations by dominated convergence on the finite history. Therefore

    |K|<=T*B1+3B0,
    |K'|<=T*B2+3B1+|h|*B0,
    |K''|<=T*B3+3B2+2|h|*B1+|h'|*B0.

Substituting the derived A.7 B_j yields precisely L0,L1,L2 in the
formulation. The exact kernel API uses the full scale factor and actual
conformal coordinate. Time differences may be evaluated as integrals in
y; none is replaced by eta_star times a background-coordinate difference.

## 3. Euler contribution at the same proper time

On the target plateau p=y^-4, f=1-delta/y⁴. The actual clock gives

    h=1/(eta_star*y*f^(3/4)),
    U=2*delta/(eta_star²*y⁶*f^(3/2)).

Consequently the ratios of the pure Euler contributions to the A.3
components at that proper time are, with z=delta/y⁴,

    rho_Euler/rho_ref=(1-z)^(-4),
    p_Euler/p_ref=(1+8z/5)*(1-z)^(-4),
    E_Euler/E_ref=(1+4z/3)*(1-z)^(-4).

They are increasing from one. Since delta<=1/2 and y>=2, z<=1/32.
Their derivatives are bounded by 8,12,12, respectively. Indeed their
numerators over (1-z)^5 are 4, `28/5+24z/5`, `16/3+4z`; evaluate these
upper bounds and the denominator lower bound at z=1/32. The certificate
checks the resulting strictly positive rational margins.

The mean-value inequality, `y^-12<=1/4096`, and the reference denominators
960,576,320 give the Euler error constants

    (1/491520, 1/196608, 3/327680)

in units `hbar*delta/(pi²*A⁴*eta_star⁸)`. This compares quantities at
the same original proper time t=A*eta_star²*y²/2; no perturbed-time
identification is assumed.

## 4. Full actual-stress coefficients

Split S into `S_B=-K/2+U/10` and the actual nonlinear Rmode. Let
`s_j=L_j/2+b_j/10`, so
`|S_B^(j)|<=delta*s_j/eta_star^(j+2)` for j=0,1,2. Applying the exact
A.7 stress reconstruction bound to S_B, using the same active history
and geometric caps, gives first-degree coefficients

    d_rho=[2s1+(4+delta_bar*b0+3delta_bar*b1)*s0]/64,
    d_p=[s2+6s1+(12+delta_bar*b0+3delta_bar*b1)*s0]/192,
    d_E=[s2+8s1+(16+6delta_bar*b1)*s0]/128.

Here `|integral U'*S_B|<=3*delta²*b1*s0/eta_star⁴`; replacing its
extra delta by delta_bar gives a bound valid throughout the closed
amplitude interval, not only at its endpoint. Add the Euler constants
from Section 3 to obtain exact C1.

The nonlinear Rmode terms give exactly A.7's actual-minus-Born error,
including its own U'*Rmode history. For the remaining R² trace piece,
`|L|<=6*delta²*b0²/eta_star⁴` on the active history, and
`a^-4<=1/(8*A⁴*eta_star⁴)` on target. Its density, pressure and EED
coefficients are `3*b0²/128`, `7*b0²/768`, `13*b0²/512`, respectively.
Add these to A.7's C2 coefficients. This proves the full actual/reference
estimate, with no uncomputed quantum response input.

All coefficients and rounded upper bounds are checked both with SymPy
exact arithmetic and a separate Fraction-only assembly from the pinned
A.7 report. No numerical sampling or rounded floating comparison proves
an inequality. The zero-amplitude endpoint holds by the exact identities.

## 5. Finite actual SEE residual and denominator 1922

The quantum source is epsilon times the actual stress. On the target,
the original prepared metric equals the order-reduced A.5 metric, so the
exact physical frozen-reference defects are

    D_i,frozen=c_i*z²*(2-z)/(4*t²*(1-z)²), c_i=(3,5,9).

The metric is used here; the surrogate fluid is not substituted for actual
RSET. With t=A*eta_star²*y²/2 and z=delta/y⁴,

    D_i,frozen=c_i*delta²*(2-z)/(A²*eta_star⁴*y¹²*(1-z)²).

For z<=1/32,

    (2-z)/(1-z)² <= 2*(32/31)² =2048/961.

An explicit polynomial positivity certificate sets v=1-32z and expands
`(2048/961)*(1-z)²-(2-z)` into

    1/32+(97/992)*v+(2/961)*v²,

whose coefficients are positive for v>=0. Combining y¹²>=4096 gives
`2048/(961*4096)=1/1922`, literally `1922=2*31²`. This denominator
is not a power or typographical shorthand. The equivalent expression
`[(1-z)^(-2)-1]/z` has continuous value 2 at zero.

Now `D_i,actual=D_i,frozen-epsilon*kappa*(T_i,actual-T_i,ref)` and

    epsilon*kappa*hbar/pi²=2880*delta*A²*eta_star⁴

by the original definition of d and delta. The triangle inequality gives
the stated finite residual. Both the physical and FK sign dictionaries
are checked; taking absolute values does not change this conclusion.

## 6. Physical finite example and honest denominator choices

The error/reference ratios are bounded by
`den_i*3⁸*(delta*C1_i+delta²*C2_i)`, where den_i=960,576,320.
These nonnegative polynomials increase with delta. At delta=10^-14 the
rounded rational bounds are all below 1/100; therefore the same is true
throughout 0<=delta<=10^-14. This proves positivity only of the actual
specified transported-state components, since their references are positive.

For the residual comparison use the positive **zeroth-order classical**
Einstein/radiation components `(3,1,3)/(A²*eta_star⁴*y⁴)`. The uniform
ratios therefore use multipliers 27,81,27 times delta² and the displayed
residual brackets. Their exact rational values at delta=10^-14 are all
below 10^-17. Neither comparison divides by unknown actual values.

Finally, `delta=4*epsilon*d/t_star²`; at physical epsilon=1,
`t_star>=2*10^7*sqrt(d)` implies delta<=10^-14 exactly. This is a
dimensionally explicit regime for the given preparation, not a statement
about observed cosmological scales or a singularity time.

## 7. What this does not establish

The earlier past metric is deliberately off shell; the present residual
estimate is on the target plateau only. A small residual does not supply
self-consistent initial data, existence, uniqueness, stability or a nearby
SEE solution. There is no change to the original prepared metric, no
certificate for the different second-order corrected family, and no
full two-point-function/QSEI estimate. In particular positive EED in this
one transported state is not a lower bound for all Hadamard states.
