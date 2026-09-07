# History-based photon-QSEI theorem and its exact quantifiers

## 1. The actual field/source inequality

The model, convention and finite prescription are exactly those of the
pinned A.16 Maxwell branch. In particular R_UU=3(Hdot+H^2), K=3H(0),
and an actual solution of G+Lambda g=-kappa(T_M+T_other) obeys

    R_UU=Lambda-kappa(E_M+E_other).

The A.16 all-Hadamard result applies to real compact proper samplers on
any smooth strip segment, and to H2_0 with both zero boundary traces:

    integral E_M f^2 >= -hbar/(8pi^2) integral |L_H f|^2
                       +integral E_conf f^2,
    L_H f=f''-2Hf'+(3H^2/4-3Hdot/2)f.

Its target need not be spatially homogeneous or extend as a state to global
Minkowski spacetime. Those facts are not changed here. The chosen actual
SEE source is one of these admitted states; an all-target field inequality
does not assert that every target state solves the same SEE.

Let Sigma=max(Lambda,0)+kappa max(-ell,0), with E_other>=ell. Then

    integral R_UU f^2 <= Q[f],
    Q[f]=kappa hbar/(8pi^2) integral |L_H f|^2
          -kappa integral E_conf f^2 +Sigma integral f^2.

Both source coarsenings are upwards; no sign is assigned to unspecified
matter. The A.16 E_conf is retained with its independent real beta_M. Its
absolute magnitude is bounded by hbar V/(2880pi^2) from the given proper
Hubble jets. An additional left-hand +alpha I term, if present, would add
alpha I/kappa to the effective additional source; its EED is not dropped.

## 2. Matched two-sided sampler and the past square

Write p(x)=3x^2-2x^3. Put

    u(t)=p((t+r tau)/(r tau)), -r tau<=t<=0,
    g(t)=p((tau-t)/tau), 0<=t<=tau,
    f=u on the past, f=g on the future, zero elsewhere.

The values match to 1 at 0 and both first derivatives vanish there. Values
and first derivatives vanish at the outer endpoints, so f is C1, piecewise
smooth and H2_0. No C-infinity claim is made for the piecewise polynomial;
the pinned H2 density result supplies the QSEI. It is used only under the
contradiction assumption that the smooth metric reaches tau and the past
interval actually exists, so its support lies compactly in a smooth strip.
No sampler is passed through a singular or absent endpoint.

For u on the past, direct differentiation yields

    R_UU u^2=3(u'-Hu)^2-3u'^2+(3Hu^2)'.

Its endpoint term is K, because u(-r tau)=0 and u(0)=1. Thus, writing
J[g]=integral_future(3g'^2+R_UU g^2),

    J[g]+K <=3(||g'||^2+||u'||^2)+Q[f]
               -3 integral_past(u'-Hu)^2.                    (1)

This is an exact past geometric quantity, not an assumed bound rho_min
on R_UU. The identity is the FLRW analogue of retaining, rather than
coarsening, the past term in Fewster–Kontou Scenario 2. It removes the
pointwise SEC premise; it does not remove initial geometric information.

On the past u,u'>=0 and H<=-h/tau imply

    3 integral_past[(u'-Hu)^2-u'^2]
      >= (6h/tau) integral u u' +(3h^2/tau^2) integral u^2
      = [3h+39h^2r/35]/tau.                                  (2)

Here integral u u'=1/2 and integral u^2=13r tau/35. Also
3||g'||^2=18/(5tau). H(0)<0 alone does not supply the quantitative duration
r tau: continuity gives some local interval, but its length must still be
verified. The code's history counterexample keeps H(0) fixed while changing
the weighted history, preventing an instantaneous-K substitution.

## 3. Exact rational cost of the past and relative future bounds

All norms here are proper-time L2 norms on the indicated sampler piece.
The unit-cubic moments are

    integral p^2=13/35, integral p'^2=6/5, integral p''^2=12,
    integral (p'/s)^2=12, integral (p/s^2)^2=13/3.              (3)

The last two expressions extend polynomially to s=0. This is why a
quadratic endpoint zero, not a linear cutoff, tolerates the relative caps.

On the past, |H^(j)|<=d_j/tau^(j+1). Set q_d=3d0^2/4+3d1/2. Triangle
inequality gives

    tau^(3/2)||L_H u||
      <= sqrt(12)r^(-3/2)+2d0 sqrt(6/5)r^(-1/2)
         +q_d sqrt(13/35)r^(1/2)
      <= [7/2+(11/5)d0 r+(2/3)q_d r^2]/r^(3/2).             (4)

The squared rational constants are strictly larger than the radicals:
(7/2)^2>12, (11/10)^2>6/5, (2/3)^2>13/35. The reference loss contributes
13r V(d)/(12600 tau^3) in units kappa hbar/(8pi^2), because 8/2880=1/360.

On the future set s=(tau-t)/tau. The relative bounds and (3) give

    ||g''||=sqrt(12)/tau^(3/2),
    ||H g'||<=c0 sqrt(12)/tau^(3/2),
    ||H^2 g||<=c0^2 sqrt(13/3)/tau^(3/2),
    ||Hdot g||<=c1 sqrt(13/3)/tau^(3/2).

Consequently

    tau^(3/2)||L_H g||
      <=(7/2)(1+2c0)+(25/12)(3c0^2/4+3c1/2),                (5)

where (25/12)^2>13/3. Each term in the finite-beta reference EED has
homogeneous weight four. Hence its magnitude is bounded by
hbar V(c)/[2880pi^2(tau-t)^4], and the reference loss is
13V(c)/(1080tau^3) in the same units. The cubic's quadratic zero makes
that weighted integral finite; a linear zero would fail this calculation.

Equations (4),(5) establish the C_P,C_F in the formulation, with beta_M
visible through V. The source norm is integral f^2=13tau(1+r)/35. Thus

    tau Q[f] <=delta(C_P+C_F)+13(1+r)sigma/35.                (6)

All past and future caps are hypotheses. In particular future relative
bounds do not assert a uniform finite H or Hdot up to tau and are not the
absolute cap Hmax*tau<=1 excluded in A.16. They permit curvature growth
like successive inverse powers of tau-t. The estimate is applied on a
hypothetically smooth segment reaching tau, where the sampler and QSEI
are unproblematic; an improper integral through a physical singularity
is not used. Neither this argument nor the field QSEI derives these caps
from initial data or controls a quantum-gravity endpoint.

## 4. Global FLRW implication without a compact-Cauchy import

Combining (1),(2),(6) gives J[g]+K<=-M/tau with M as stated. Suppose a
comoving normal reaches tau. Because the global metric is smooth and a>0,
this normal is nonfocal: its transverse Jacobi scale is a(t)/a(0)>0. More
directly the same square identity on the future gives

    J[g]+K=3 integral_0^tau (g'-Hg)^2 dt >0.                 (7)

Strict positivity follows without a pointwise Ricci sign: equality would
force g'=Hg almost everywhere and hence g(t)=a(t)/a(0), which cannot obey
g(tau)=0 while a(tau)>0. Equation (7) contradicts M>=0.

Thus no comoving normal reaches tau. In the stipulated global product
spacetime, this means the upper endpoint t_+ of I_t is <=tau. The comoving
normal t in [0,t_+) is inextendible within M and has finite proper length.
For any future timelike curve from t=0, ds<=dt, so its proper length is
at most t_+<=tau as well. No compact-Cauchy theorem is invoked verbatim:
the spatial R^3 is noncompact and the global FLRW clock provides this
direct conclusion.

The past interval is not conditional; it must really exist. The future
cap hypothesis is only about a normal segment reaching tau. If a curve
ends earlier, the conclusion already holds and no geometric continuation
across its endpoint is manufactured. One must not strengthen the certified
claim to say all shorter segments satisfy these relative caps, or weaken
it by assuming an unverified future cap from a short initial patch.

Incompleteness is a statement about M. The proof does not establish that M
admits no larger extension or that a coordinate-patch endpoint is a
fundamental singularity. If the intended physical spacetime is a larger
extension, the theorem's global hypotheses must be checked there. Time
reversal replaces H by -H(-t) and reverses the normal; the squared field
and geometric functionals give the corresponding past conclusion.

## 5. Rational macroscopic calibration and source budget

For r=1/100,h=3/2,d=c=(2,4,16,96), the worst |beta_M|<=1 value gives
V=16704 and the exact C_P,C_F in the formulation. Their sum is strictly
less than 13,000,000. At delta=10^-8 the past gain is 63351/14000, the
future gradient is 18/5, and the remaining zero-source margin is

    2398243151869/3000000000000 >3/4.

The source weight is 1313/3500. Thus sigma<=1 leaves more than
3/4-1313/3500=328/875>1/3. These calculations are exact; no floating-point
search or approximate inverse is used. Coefficients and signs are checked
by an independent Fraction engine and a separate covariant audit.

For kappa=8piG the quantum parameter is ell_P^2/(pi tau^2). The condition
is therefore compatible with very long macroscopic and cosmological tau.
This is a scale statement within the conditional theorem, not a verified
observational application, a thermal-radiation witness, or an assertion of
EFT validity at arbitrarily large curvature. No fixed beta is inherited
from the old scalar field; the generic theorem keeps it explicit.

## 6. Smooth future-complete geometric countercontrol

This control prevents the geometric hypotheses from being mistaken for
an independent incompleteness theorem. In x=t/tau define p5(x)=10x^3-
15x^4+6x^5 on [0,1], extend p5 by0 before0 and1 after1, and put

    H_raw(t)=-(31/(20tau))[1-p5(t/tau)].

It is C2 with bounded weak third derivative, that is W^{3,infinity}.
Convolve in t with any positive, normalized C-infinity kernel supported
in [-tau/1000,tau/1000]. The resulting H is smooth and its derivatives
through three are convolutions of the bounded weak derivatives. Positive
unit-L1 convolution cannot increase their supremum norms.

On [0,1], p5'=30x^2(1-x)^2<=15/8. Writing w=(2x-1)^2 gives
p5''^2=225w(1-w)^2 and

    4/27-w(1-w)^2=(w-1/3)^2(4/3-w)>=0.

Thus |p5''|<=sqrt(100/3)<6. Finally p5'''=60-360x+360x^2 lies between
-30 and60. The resulting global dimensionless caps for H are

    (31/20,93/32,93/10,93) < (2,4,16,96).

These imply all the relative future caps on 0<=t<tau. For t<=0 the
convolution samples raw future x at most1/1000; since 0<=p5(x)<=10x^3,

    -tau H(t)>=(31/20)(1-10^-8)>3/2.

Hence the entire required past history also holds. H' is nonnegative,
so R_UU=3(H'+H^2)>0 on that short past: the pointwise SEC is violated,
not hidden in the history assumption.

Define a(t)=a(0)exp(integral_0^t H). For t>=0, a is positive, decreases
by at most exp[-(31/20)(1001/1000)], and becomes constant by t=(1001/1000)tau.
Every timelike geodesic has conserved spatial momentum and proper-time
element dt/sqrt(1+P^2/a^2), bounded below by a positive multiple of dt
on the future. Null affine length is proportional to integral a dt.
Both diverge as t tends to infinity, proving future causal completeness.

This is NOT an allowed small-delta Maxwell SEE solution. For example,
in its constant past the actual conformal Maxwell vacuum has I=0 and
kappa tau^2 E_conf=-delta(31/60)(31/20)^4. At Lambda=0 its required
negative additional EED would have dimensionless magnitude

    3(31/20)^2-delta(31/60)(31/20)^4 >1 at delta=10^-8,

already outside the sourced calibration budget. More generally, the
theorem excludes any permitted actual SEE/Maxwell state matching this
complete geometry with those small quantum/source constants. The control
demonstrates that the field/source inequality matters; it is not a false
state witness or a second physical cosmology. It also avoids calling a
contracting de Sitter flat patch future causally complete.
