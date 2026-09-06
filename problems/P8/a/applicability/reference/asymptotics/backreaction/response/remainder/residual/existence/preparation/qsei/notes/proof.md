# Exact scattering, two-frequency control and the absolute QSEI

Every geometry and state in this proof is the actual A.11 solution on
the stated slab and its unchanged original past. Smoothness is supplied
by A.11. Constants below use only quantitative C1 potential information.

## 1. Positive-type input and exact rescaling

For this massless minimal field, a difference of EEDs is a difference
of squares of the proper-time field derivative. State-independent
renormalization terms cancel in that difference. The differentiated
target and reference two-point functions are positive type, have the
same commutator, and their difference is smooth and symmetric. Their
timelike pullbacks exist by the Hadamard wavefront condition.

Apply the positive-type diagonal Fourier argument, with its one-sided
factor 1/pi, to the pulled-back kernel. It holds for every Hadamard
target state, not just quasifree or homogeneous ones. This is the
specialization of FK Lemma 4 audited in immutable A.9; its geometric
constants are not imported. Coordinate rescaling and multiplication by
real smooth weights preserve positivity. If psi(s)=f(T0*s) and

    F(x)=a(x)^(-3/2) psi(s(x)), D=partial_x-h(x),

the exact difference inequality is

    integral f^2(E_omega-E_reference) d_tau >= -Q[f],
    Q[f]=hbar/(4 pi^3 T0^3) I[F],
    I[F]=integral_0^infinity k dk integral_0^infinity d_alpha
              |integral e^(-i alpha x) F(x) D w_k(x) dx|^2.

The modes have original past value `w_k=e^(-ikx)` and Wronskian
normalization corresponding to `w_k=sqrt(2 k_comoving) v_k` before
the dimensionless momentum change. The reference radial measure is
`hbar k dk/(4 pi^2)`. Constant phase changes from shifting the origin
of x cancel in the spectral squared norm. No inverse proper-clock
expansion occurs; `ds=a dx` is used exactly. The estimates below prove
finiteness of the reference positive-frequency integral directly.

All unlabelled time L2 norms below are in dx; the potential suprema
hold on its entire history, while h, h', h'' suprema need hold only
on the sampler interval. The Fourier convention is
`hat v(omega)=integral exp(-i omega x) v(x) dx`, with full Parseval
`integral_R |hat v|^2=2 pi integral_R |v|^2`.

## 2. Exact scattering split with C1 potential bounds

Use history duration T=3, derivative cap M=10^-5 and
B=MT=3*10^-5. Since u(0)=0, `|u|<=B` on this whole history.
The actual modes satisfy

    b_k=e^(ikx) w_k,
    b_k'=-integral_0^x exp(2ik(x-r)) u(r)b_k(r) dr,
    b_k-1=-integral_0^x [(exp(2ik(x-r))-1)/(2ik)] u(r)b_k(r) dr.

The Volterra modulus estimate gives
`|b_k|<=cosh(sqrt(B)T)<=exp(BT^2/2)<2`, for every k>0.
Indeed `(2n)!>=2^n n!` proves the exponential estimate termwise.
This is valid without assuming k^2+u positive.

Write q_k=u b_k, and define the actual complex histories

    A_k(x)=integral_0^x q_k(r) dr,
    B_k(x)=integral_0^x exp(-2ikr) q_k(r) dr.

The repeated letter B in B_k denotes a function; the un-subscripted
B remains the real potential cap. Then the exact differentiated mode is

    D w_k=e^(-ikx)(-ik-h)+e^(-ikx) c_f A_k+e^(ikx)c_b B_k,
    c_f=-1/2-h/(2ik), c_b=-1/2+h/(2ik).

One integration by parts, with **q_k(0)=0**, gives

    B_k=-e^(-2ikx)q_k(x)/(2ik)
           +(1/(2ik)) integral_0^x e^(-2ikr) q_k'(r) dr.

Had q_k(0) been nonzero, an additional q_k(0)/(2ik) term would be
present. This endpoint is explicitly retained in the symbolic audit.

The elementary bounds are

    |b_k'|<=2BT, |b_k-1|<=BT^2,
    |q_k|<=Q0=2B, |q_k'|<=Q1=2M+2TB^2,
    |A_k|<=A0=2BT, |A_k'|<=A1=Q0, |A_k''|<=A2=Q1.

For k>=1 the two Volterra identities and the same single integration
by parts sharpen two of these:

    |b_k-1|<=l0/k, l0=2BT,
    |b_k'|<=l1/k, l1=B+TM+T^2 B^2.

Here `|b_k'|=|B_k|` and the bound on q_k' just proved is sufficient.
There is no derivative of u' in this argument. Consequently

    q_k'=u'+r_k, |r_k|<=R1/k, R1=M l0+B l1.

This leaves precisely four ultraviolet error pieces:

    forward:       e^(-ikx)c_f A_k,
    local:        -e^(-ikx)c_b q_k/(2ik),
    real history:  e^(ikx)c_b/(2ik) integral_0^x e^(-2ikr) u'(r) dr,
    remainder:     e^(ikx)c_b/(2ik) integral_0^x e^(-2ikr) r_k(r) dr.

The identity is exact, not a Born approximation. The r_k remainder
still contains all nonlinear state evolution.

On the new sampling interval, h<=H0=1/2. The identities
`h'=-u-h^2`, `h''=-u'-2hh'` and A.11's local u and global u' caps
give `|h'|,|h''|<=H1=H2=1/3`. Thus the forward amplitude has two
derivative bounds

    C0=(1+H0)A0/2,
    C1=[(1+H0)A1+H1 A0]/2,
    C2=[(1+H0)A2+2H1 A1+H2 A0]/2.

The local piece is `e^(-ikx) P_k/k`, where

    |P_k|<=P0=(1+H0)Q0/4,
    |P_k'|<=P1=(1+H0)Q1/4+H1 Q0/4.

For both history pieces `|c_b/(2ik)|<=A/k`, where A=(1+H0)/4=3/8.
These yield exact constants

    (C0,C1,C2)=(27/200000,3/40000,1300081/20000000000),
    (P0,P1)=(9/400000,500081/40000000000),
    R1=3600243/10^15.

## 3. Integrating both frequencies without higher jets

The conformal sampler envelope has length at most ell=L/2. Two
sampler integrations by parts in the forward piece give

    I_forward <= (ell/3) [C0||F''||+2C1||F'||+C2||F||]^2,

because the exact double moment is
`integral_1^infinity k dk integral_0^infinity(alpha+k)^(-4) d_alpha=1/3`.
One sampler integration by parts in the local piece suffices:

    I_local <= ell [P0||F'||+P1||F||]^2,

using the double moment with weight `1/[k^2(alpha+k)^2]`, which is 1.
Both integrations have zero sampler boundary terms.

For the real history put
`J_k(x)=integral_0^x exp(-2ikr) u'(r) dr`. First use **full** Parseval
in the sampling frequency. The product F times the history is generally
complex, so equality with half of that norm is not justified. Tonelli
for the resulting nonnegative expressions gives

    I_history <= 2 pi A^2 integral dx |F(x)|^2
                              [integral_1^infinity |J_k(x)|^2 dk/k].

The following uniform estimate removes the inner integral. The truncated real function
`u' 1_[0,x]` is L2 even though it has an endpoint jump. Parseval in r,
the reality symmetry, and the frequency change omega=2k imply

    integral_0^infinity |J_k(x)|^2 dk
        =(pi/2) integral_0^x |u'(r)|^2 dr <=(pi/2)T M^2.

Since 1/k<=1 on k>=1, this proves

    I_history <= pi^2 A^2 T M^2 ||F||^2.

This second use of Parseval is the essential improvement over demanding
a pointwise k^-1 bound on J_k from u''. The first Parseval was full;
the second uses half-frequency symmetry only for the explicitly real
history. Neither is replaced by a frequency cutoff calculation.

The nonlinear backward remainder is pointwise <=A T R1/k^2, so full
sampling Parseval and `integral_1^infinity k^-3 dk=1/2` give

    I_remainder <= pi A^2 T^2 R1^2 ||F||^2.

For 0<k<1 no oscillatory integration by parts is needed. The total
mode error obeys `|e_k|<=E_IR=(1+H0)BT^2+2BT=117/200000`.
Full Parseval and `integral_0^1 k dk=1/2` give

    I_IR <= pi E_IR^2 ||F||^2.

All exchanges here are justified by nonnegative Tonelli or by the
finite norm bounds just displayed. The spectral triangle inequality
combines the four ultraviolet pieces and the disjoint infrared error.

## 4. Auxiliary norm, proper sampler and rational coefficient

The auxiliary differentiated mode `e^(-ikx)(-ik-h)` is a comparison
kernel only; it is not claimed to define a state of the perturbed
equation. Its exact spectral norm, for the actual real h(x), is

    I0=(pi/4) integral [F''^2+2((hF)')^2+(8/3)F''(hF)'] dx,
    sqrt(4I0/pi)<=||F''||+(3/2)||(hF)'||.

The integrand has the positive-square decomposition recorded in
sampling.py. This identity is rechecked; no old numerical coefficient
is used. Combining all pieces gives

    sqrt(4I/pi) <= ||F''||+(3/2)||(hF)'||
      +(2/sqrt(pi))*sqrt(ell/3)*(C0||F''||+2C1||F'||+C2||F||)
      +(2/sqrt(pi))*sqrt(ell)*(P0||F'||+P1||F||)
      +2sqrt(pi) A M sqrt(T)||F||+2 A T R1||F||+2E_IR||F||.

Let D2=||psi_ss|| in ds. Since ds=a dx, the proper envelope has
length <=3ell. Two Dirichlet Poincare estimates with pi>3 give
`||psi_s||<=ell D2`, `||psi||<=ell^2 D2`. On the actual metric
`H=h/a<=1/4` and `|H_s|=|u+2h^2|/a^2<=1/7`.
The exact proper-clock identities, including the square-root measure
factor, then give

    ||F||<=f0 D2,        f0=ell^2/4,
    ||F'||<=f1 D2,       f1=(ell+3ell^2/8)/2,
    ||F''||<=f2 D2,      f2=1+ell/2+117ell^2/448,
    ||(hF)'||<=b1 D2,    b1=ell/4+39ell^2/224.

Explicitly `F''/sqrt(a)=psi_ss-2H psi_s+(3H^2/4-3H_s/2)psi`
and `(hF)'/sqrt(a)=H psi_s+(H_s-H^2/2)psi`. These are checked in
the independent clock algebra. They are not background-clock formulas.

Use only upward rational replacements: `sqrt(ell),sqrt(ell/3)<=10^-5`,
`2/sqrt(pi)<2`, `2sqrt(pi)<4` and `sqrt(T)<2`. The resulting exact root is

    6400000017560000000992361942917377673
    /6400000000000000000000000000000000000 <1+10^-8.

Its square is strictly below 2. The exact arithmetic is independently
reassembled with Fraction, and the individual root contributions are
reported rather than hidden in this final quotient. Finally
`||f_tau_tau||^2_d_tau=T0^-3 D2^2`, so

    Q[f] <= 2 hbar/(16 pi^2) ||f_tau_tau||^2_d_tau.

## 5. New positive credit from the actual unforced SEE

On the free half the preparation source is zero. The full actual SEE
with the fixed ordinary radiation therefore gives

    E_reference= hbar/(pi^2 A^4 eta_star^8)
                    *[a^2(u+h^2)-1]/(960 delta a^4).

This is not a statement about arbitrary target states, and is not an
old-metric RSET value inserted into the new geometry. The pointwise
equality is valid because A.11 already proved the smooth full equation.

The old plateau at the same x has numerator `3z/(1-z)`, z=delta/y^4,
which is >=delta/27 because y<=3. A.11's full contraction sharpens its
ball radius to `||X*||<rstar=72*10^-9`. Its exact pair bounds imply

    |Delta[a^2(u+h^2)-1]|
       <=[9L+(9/2)L^2+(3B_local+3/4)L^3] rstar <delta/54,
    B_local=2*10^-12.

Hence the new numerator is >=delta/54. Since a^4<=81, the dimensionless
reference EED is >=1/4199040. Dropping this newly proved positive credit
from the difference inequality yields the claimed absolute bound.

## 6. States, density and limitations

Every estimate is uniform along comoving lines by reference homogeneity;
no homogeneity requirement is imposed on the target state. For each
fixed target Hadamard state the renormalized EED is smooth on compact
segments. Smooth compact samplers are dense in H2_0 there, the spectral
norm is H2-bounded, and the proper/conformal changes are smooth and
nonsingular. Thus the inequality extends to the stated H2_0 domain.

There is no new field-independent or massive/nonminimal QSEI. The
absence of u''/u''' bounds here comes from a different two-frequency
argument, not from interpreting qualitative smoothness as a uniform
higher-jet estimate. The short sampling domain and the fixed absolute
renormalization prescription remain essential. The incomplete focusing
application is quantified separately, rather than silently promoted to
a Penrose--Hawking conclusion.
