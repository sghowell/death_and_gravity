# Fresh signed-reference, all-sampler QSEI and actual index estimate

All geometry below is the actual A.14 smooth solution. The source/state,
initial data, renormalization scheme and the original preparation cutoff
are fixed as in [FORMULATION.md](../FORMULATION.md). Primes on a,h,u are
dimensionless conformal derivatives; proper dimensionless time satisfies
ds=a dx. Set ell=L=10^-6, a safe envelope for the whole source-free domain
I=(x0+L0/2,x0+L), where L0=10^-10. Its exact conformal width is L-L0/2.

## 1. New geometry supplies the general C1 mode hypotheses

A.14 proves a smooth positive metric and the original vacuum transported
by its Klein--Gordon evolution. In particular it proves an actual
Hadamard reference, not a finite-frequency mode surrogate. Its full
original history still has duration T<=3 and |u'|<=M=10^-5. Since
u=0 in the original free past, |u|<=B=MT=3*10^-5 on that history.
On the new sampling domain it proves

    2<=a<=3, 1/3<=h<=1/2, |u|<=3*10^-12.

The exact identities h'=-u-h^2 and h''=-u'-2hh' therefore give
|h'|,|h''|<1/3, with fresh strict margins. These, and only these
quantitative potential/Hubble hypotheses, are needed for the general
C1 scattering estimates proved in A.12 notes/proof.md sections 2–3.
No unknown quantitative higher jet is supplied by qualitative smoothness.

For clarity, that lemma is an infinite-frequency estimate: for the
actual modes w_k=e^(-ikx)b_k with the original free-past normalization,
q_k=u b_k and D=partial_x-h, the exact forward/backward split is

    D w_k=e^(-ikx)(-ik-h)
       +e^(-ikx)c_f integral_0^x q_k(r)dr
       +e^(ikx)c_b integral_0^x e^(-2ikr)q_k(r)dr,
    c_f=-1/2-h/(2ik), c_b=-1/2+h/(2ik).

One backward integration by parts retains the endpoint q_k(0)=0;
q_k'=u'+r_k, |r_k|<=R1/k for k>=1. The Volterra majorant proves
|b_k|<2 for all k>0 from exp(BT^2/2)<2. The derivative of u' is not
needed. The backward real history is controlled by momentum Parseval
for u' 1_[0,x], while the sampler Parseval is **full**, because its
product with the history is complex. This distinguishes the two
frequency arguments and retains the entire nonlinear state response.

Applying the general majorants with the **new A.14 inputs** yields

    (C0,C1,C2)=(27/200000,3/40000,1300081/20000000000),
    (P0,P1)=(9/400000,500081/40000000000),
    Aback=3/8, R1=3600243/10^15, EIR=117/200000.

These happen to equal A.12's mode constants because its history caps
still hold; neither its sampler length nor its reference credit does.
Both the inputs and the new local h'/h'' margins are independently
recomputed here, and the pinned general symbolic identities are replayed.

## 2. Positive type and the new proper sampler interval

Let psi(s)=f(T0 s), F(x)=a(x)^(-3/2)psi(s(x)). For every Hadamard target
state the positive-type diagonal Fourier argument, with the one-sided
factor 1/pi, gives the exact difference inequality

    integral f^2(E_target-E_reference)d_tau >= -Q[f],
    Q[f]=hbar/(4 pi^3 T0^3) I[F],
    I[F]=integral_0^infinity k dk integral_0^infinity d_alpha
         |integral e^(-i alpha x) F(x)D w_k(x)dx|^2.

This uses smoothness and equal commutator of the Hadamard difference,
positive type of the differentiated two-point functions and the
timelike pullback wavefront condition. It is the general input audited
in A.9/A.12, not a new restriction to quasifree target states. For this
massless minimal field the EED difference is the proper-derivative
square difference; state-independent local terms cancel.

The general two-frequency estimates, now on a conformal envelope ell,
are

    I_forward <= (ell/3)(C0||F''||+2C1||F'||+C2||F||)^2,
    I_local <= ell(P0||F'||+P1||F||)^2,
    I_history <= pi^2 Aback^2 T M^2 ||F||^2,
    I_remainder <= pi Aback^2 T^2 R1^2 ||F||^2,
    I_IR <= pi EIR^2 ||F||^2.

The first two integrate the sampler and momentum variables, with exact
moments 1/3 and 1. The real-history identity is
`integral_0^infinity |J_k(x)|^2 dk=(pi/2)integral_0^x |u'|^2`,
where J_k is the Fourier transform of the truncated real u' at 2k.
The separate full sampling Parseval contributes 2pi, not pi.
Nonnegative Tonelli and these finite norms justify the exchanges.

The auxiliary differentiated mode e^(-ikx)(-ik-h) is only a comparison
kernel. Its exact norm satisfies

    sqrt(4 I0/pi)<=||F''||+(3/2)||(hF)'||.

The actual proper clock has H=h/a<=1/4 and
|H_s|=|u+2h^2|/a^2<1/7. The proper sampler envelope is at most 3ell;
Dirichlet Poincare, applied twice with pi>3, gives
`||psi_s||<=ell||psi_ss||`, `||psi||<=ell^2||psi_ss||`.
The exact measure-adjusted clock identities are

    F/sqrt(a)=psi/a^2,
    F'/sqrt(a)=(psi_s-(3/2)H psi)/a,
    F''/sqrt(a)=psi_ss-2H psi_s+(3H^2/4-3H_s/2)psi,
    (hF)'/sqrt(a)=H psi_s+(H_s-H^2/2)psi.

Consequently, in units D2=||psi_ss||_(ds), the four relevant norms are
bounded by

    f0=ell^2/4,
    f1=(ell+3ell^2/8)/2,
    f2=1+ell/2+117ell^2/448,
    b1=ell/4+39ell^2/224.

Use the new exact sqrt(ell)=1/1000, sqrt(ell/3)<=sqrt(ell),
2/sqrt(pi)<2, 2sqrt(pi)<4 and sqrt(T)<2. The total spectral root is

    R=f2+(3/2)b1
       +(2/1000)(C0 f2+2C1 f1+C2 f0)
       +(2/1000)(P0 f1+P1 f0)
       +8 Aback M f0+2 Aback T R1 f0+2 EIR f0
      <1+2*10^-6.

The inequality is checked as an exact rational, including every error
term. It is about 1.000001145; that decimal is not certificate input.
Since ||f_proper_second||^2=T0^-3 D2^2,

    Q[f]<=hbar/(16 pi^2) R^2 ||f_proper_second||^2.

## 3. Actual extended reference: signed, not assumed positive

The preparation source is identically zero on I. A.14 already proves
the pointwise full SEE, with its fixed ordinary radiation. Subtracting
that radiation from the Einstein effective energy gives

    E_reference=hbar/(pi^2 T0^4) N/(960 delta a^4),
    N=a^2(u+h^2)-1.

This equality is only for the actual reference SEE source, not arbitrary
target states. During the preparation, its external EED would also
have to be subtracted; the support restriction is essential.

Compare with the unchanged prepared baseline at the same x. Its
numerator is Nbar=3delta/(y^4-delta)>0. For the actual A.14 fixed point,

    ||X*||_sigma<rstar=9*10^-8, sigma=2*10^6,
    exp(sigma L)<E=9, alpha=1/sigma.

The same-data damped metric estimates in A.14 give the weighted gains
alpha, alpha^2, 18alpha^3 for Delta u, Delta h, Delta(a^2).
The exact algebraic decomposition

    N-Nbar=a^2[Delta u+(h+hbar)Delta h]
                  +Delta(a^2)(ubar+hbar^2)

therefore implies

    |N-Nbar|<=loss
       =E rstar[9alpha+9alpha^2+18(3*10^-12+1/4)alpha^3].

Both actual and baseline geometries lie in the same A.14 enclosure.
Dropping only the positive Nbar, and using a^4>=16, proves

    E_reference >= -hbar/(pi^2 T0^4) loss/(960delta*16)
                 >= -hbar/(40pi^2 T0^4).

The last coarsening is strict and recomputed with Fraction. The
intermediate negative magnitude is about 0.023731, below 1/40.
The loss is larger than the old positive-numerator guarantee; no
positive reference credit is inferred on the extended interval.
Neither the large weighted radius nor a naive unweighted length
replacement is substituted for the sharp fixed-point metric gains.

Poincare in the **new** proper envelope now gives

    integral f^2 d_tau <= T0^4 ell^4 integral |f_proper_second|^2 d_tau.

Thus the possible negative reference part contributes only

    [hbar/(16pi^2)]*(2/5)ell^4 ||f_proper_second||^2.

Exact recomputation proves `R^2+(2/5)ell^4<2`. Combining this signed
reference term with the difference inequality proves the absolute
coefficient2 stated in the formulation. The same rounded coefficient
as A.12 is an output of this new calculation, not an inherited assumption.

For each target Hadamard state the renormalized EED is smooth on compact
segments. Smooth compact real samplers are dense in H2_0 there; the
spectral norm is H2 bounded and the smooth proper/conformal change is
nonsingular. Passing to the limit proves the asserted H2_0 extension.
The homogeneous reference makes the constants uniform on comoving lines
without requiring homogeneous target states.

## 4. Fresh focusing test uses actual curvature

Let tau_hat=tau/T0. Every segment in I has tau_hat<=D=3ell=3*10^-6.
In the fixed FK convention,

    R_UU T0^2=-3(u+h^2)/a^2,
    |R_UU| T0^2<=3(3*10^-12+1/4)/4<1/4,
    |K|T0=3h/a<=3/4.

For any index trial g(0)=1, g(tau)=0, without assuming 0<=g<=1,
Cauchy--Schwarz from the final endpoint gives

    integral g^2 <= (tau^2/2) integral g'^2,
    integral g'^2 >=1/tau.

In either orientation the actual index form therefore satisfies

    J[g]T0 >=3/tau_hat-tau_hat/8 >=3/D-D/8
             =10^6-3/(8*10^6)>3/4>=-K T0.

The positive factor needed when inserting the second Cauchy--Schwarz
bound is explicitly checked. Hence the sufficient comoving index
trigger cannot hold on these segments even with exact curvature instead
of a QSEI estimate. The normal volume ratio a^3/a_initial^3>=8/27
supplies an independent absence-of-zero check within this slab.

The ordinary radiation has nonnegative EED. The absolute coefficient2
and actual SEE yield the sufficient geometric Q2=360delta T0^2, Q0=0.
The available-duration ratio is

    Q2/tau_available^2 >=360delta/(3ell)^2=2/5.

That is a comparison for this chosen sufficient coefficient, not an
optimal-field-constant lower bound. The reduction from 160000000 in
A.12 reflects the longer domain; it does not repair the actual index
obstruction or provide a cosmological duration. Choosing a closer-to-one
QSEI coefficient would not alter this actual-curvature conclusion.

Other hypersurfaces, boosts/null geodesics, realistic fields and unknown
long continuations are not excluded. An artificial edge of a proved
local interval is not a singularity or a maximal-development boundary.
No completeness/incompleteness conclusion for the spacetime is drawn.
