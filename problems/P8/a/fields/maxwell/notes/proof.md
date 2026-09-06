# Local Maxwell conformal QSEI, anomaly and conditional source matching

All clocks in this note are exact. The scalar field and preparation chain
are not reused as Maxwell physics. The immutable A.12 clock identity and
A.1 geometric implication are replayed only for the uses stated below.

## 1. Local physical algebra and conformal comparison

On I_eta times R^3, a smooth positive a preserves causal cones. In four
dimensions the Hodge star on two-forms is conformally invariant, so the
source-free physical Maxwell equations dF=d*F=0 and the conformally matched
commutator give the usual field-strength algebra identification with the
flat strip. There is no boundary or nontrivial flux/cohomology sector here.
The two-point tensor distributions transform by the corresponding smooth
weights; their Hadamard wavefront condition and positivity are preserved.

In orthonormal components T_UU=(sum E_i^2+sum B_i^2)/2. Its point-split
kernel is a sum of positive-type field-strength kernels. For two Hadamard
states the renormalized stress difference is the classical point-split
difference, so on the comoving curve

    (E_omega-E_conf)_g=(rho_omega-rho_vac)_flat/a^4.

This uses the state independence of the Maxwell trace anomaly: the trace
of a stress difference vanishes even though each renormalized trace need
not. The homogeneous factor in the covariant stress is a^-2 and the two
proper-unit velocities contribute a further a^-2.

An arbitrary strip state need NOT be extended globally to apply the flat
estimate. Restrict the Minkowski vacuum reference kernel to the strip.
The target-minus-reference point-split energy kernel is smooth symmetric
by Hadamard form and the common local commutator. Multiplication by a
real F(eta) times F(eta') with compact support in the strip gives a smooth
compact kernel that extends by zero to R^2. The diagonal Fourier identity
and positivity used in Fewster–Pfenning Theorem V.1 then require only
smearings supported in that strip. The discarded target contribution is
nonnegative there. No assertion of a global Minkowski target state is used.
The paper's compact-Cauchy theorem is not invoked verbatim on R^3.

The explicit restricted flat-vacuum kernel evaluates the bound with the
same Fourier integral as in FP97–102:

    Q_flat[F]=4 hbar/(2pi)^3 integral_0^infinity d_alpha
                         integral_0^infinity dk k^3 |Fhat(alpha+k)|^2
             =hbar/(8pi^2) integral |F''|^2 d_eta.

Here Fhat(v)=integral exp(-iv eta)F(eta)d_eta. Changing to v=alpha+k
gives integral_0^v k^3 dk=v^4/4. Reality supplies half of the full
Parseval integral, namely integral_0^infinity v^4|Fhat|^2=pi||F''||^2.
These are two Maxwell polarizations; the scalar half coefficient is
not imported. This is an exact established QEI, not an optimality claim.

## 2. Exact proper clock and H2 functional

Let f=f(t) be compactly supported inside the proper image of the strip,
dt=a d_eta, and F(eta)=a(eta)^(-3/2)f(t(eta)). The averaged stress
difference is integral F^2 rho_flat d_eta. Direct differentiation gives

    F''_eta/sqrt(a)=f''_t-2H f'_t+(3H^2/4-3Hdot/2)f=L_H f,
    integral |F''_eta|^2 d_eta=integral |L_H f|^2 dt.

This measure-adjusted identity is independently checked and agrees with
the general proper-clock identity pinned in A.12. No scalar mode equation,
short-interval scalar coefficient or inverse-clock approximation enters.
Thus the difference inequality in the formulation follows. Adding the
explicit conformal-reference EED in section 3 makes it absolute.

For later independent controls set p=-2H, q=3H^2/4-3Hdot/2. Compact
support (or both H2_0 boundary traces) eliminates the boundary primitive

    p f'^2+2q f f'+(p q-q')f^2.

Consequently integral (L_H f)^2 equals integral [f''^2+A1 f'^2+A0 f^2],
where

    A1=p^2-p'-2q=(5/2)H^2+5Hdot,
    A0=q^2+q''-(pq)'
      =9H^4/16+9H^2 Hdot/4+3Hdot^2/4
        -3H Hddot/2-3Hthird/2.

The difference of these densities is exactly the derivative of the
displayed primitive. Omitting the endpoint requirement is a negative
control; the identity need not hold for f=1 on a finite interval.

For each fixed target state, its renormalized EED is smooth on compact
segments. Multiplication by smooth metric coefficients and the smooth
nonsingular clock map are bounded between the relevant compact H2 spaces.
Smooth compact real samplers are dense in H2_0; the spectral bound above
is H2 continuous. Both sides therefore pass to the limit. This proves
the all-Hadamard, all-stated-H2_0 result without restricting higher point
functions, homogeneity or the one-point function of the target state.

## 3. Source/sign map and the complete finite-beta reference

Use the fixed FK convention (+---), with physical positive T_00 and

    R_00=3(Hdot+H^2), R=6(Hdot+2H^2), G_00=-3H^2.

In the Lorentzian HH normalization the metric is -+++, Ric_HH=-Ric_FK
and R_HH=R_FK; the physical covariant stress is the same tensor under
this convention change. The sign is also fixed by HH's positive physical
Casimir energy and negative mixed T^0_0 in equations 25,27,28, not chosen
from conservation alone. Define in FK

    H3_ab=R_a^c R_cb-(2/3)R R_ab
          -(1/2)g_ab R_cd R^cd+(1/4)g_ab R^2.

HH's H3_cov=-H3_FK under the simultaneous metric/Ricci convention change.
Their equation 23 is -2 a4 H3_HH/(16pi^2). For Maxwell a4=31/180 this
therefore gives +62 H3_FK/(2880pi^2), restoring hbar. The same coefficient
is checked by the Maxwell b=-62 and H3 term of Markowicz et al.; their
R=-6a''/a^3 is opposite FK, but this tensor is quadratic in curvature.

In spatially flat FLRW direct contractions give

    H3_00=3H^4, trace H3=12H^2(H^2+Hdot)=Euler4/2.

The conformal-vacuum stress follows from this tensor law and its actual
flat-vacuum stress zero. Trace plus conservation alone would leave an
undetermined C/a^4 in the density; it is the conformal state prescription
that sets that independent state term to zero. Arbitrary target stress
can have a nonzero flat counterpart, handled by the QEI, not deleted.

The type-D-free scheme is a named optional prescription. A finite R^2
counterterm shifts the stress by the independently defined conserved FK
tensor

    I_ab=2R R_ab-(1/2)g_ab R^2
         -2(g_ab Box-nabla_a nabla_b)R,
    trace I=-6 Box R.

Our sign definition is the inverse-metric variation of integral sqrt|g|R^2
in FK conventions. It is checked by direct FLRW conservation. The trace
of the new prescription is

    trace T_conf=hbar/(2880pi^2)[31 Euler4-6 beta_M Box R].

The constant beta_M is not the scalar gamma and not the c parameter of a
source written in another curvature convention. Finite Weyl-squared
variation vanishes on conformally flat metrics; the Euler variation is
topological, so after fixing the separate cosmological/Einstein terms
this one parameter describes the relevant finite curvature-squared
ambiguity. Variations are taken before restricting the metric class.

With the common factor hbar/(2880pi^2) suppressed, explicit components are

    rho=186H^4+beta_M(18Hdot^2-108H^2 Hdot-36H Hddot),
    Theta=744H^2(H^2+Hdot)
        -36beta_M(Hthird+4Hdot^2+7H Hddot+12H^2 Hdot),
    p=(rho-Theta)/3,
    E=rho-Theta/2
      =-186(H^4+2H^2 Hdot)
        +beta_M(18Hthird+90Hdot^2+90H Hddot+108H^2 Hdot).

Symbolic identities and an independent rational polynomial calculation
check every component, trace, physical sign and conservation. In de Sitter
the beta term vanishes, rho>0 and E=-rho<0: omitting the trace would change
the inequality. beta-dependent terms generally survive on other metrics.

We do NOT import the printed derivative H1 tensor of Markowicz et al.
as our I_FK. Literal conversion of its printed FLRW H1 density/pressure
under the above dictionaries gives a nonzero conservation residual
72(Hdot+H^2)(Hddot+4H Hdot). That is a specific normalization/source-use
warning, not a global claim about the paper. Its unambiguous quadratic
H3, conformal weight and Maxwell anomaly coefficient suffice here;
the finite-beta term is instead defined and checked independently.

## 4. Curvature envelopes and scheme-independent radiation control

On a compact sampler interval of proper length at most D let H_j bound
the absolute jth proper derivative of H, j=0..3. Two Dirichlet Poincare
bounds give ||f'||<=D/pi ||f''|| and ||f||<=D^2/pi^2 ||f''||. The triangle
inequality and pi>3 give ||L_H f||<=A_D||f''|| with A_D in the formulation.
The explicit stress polynomial proves E_conf>=-hbar V/(2880pi^2).
This establishes B2,B0 exactly. For beta_M=0 higher-jet terms drop; for
other schemes their numerical caps remain required physical hypotheses.

For radiation H=1/(2t), t>0, R=0 identically, I=0 and

    E_conf=31hbar/(2560pi^2 t^4),
    L_H f=f''-f'/t+15f/(16t^2).

The expanded square has coefficients -15/(8t^2),945/(256t^4). The
reference contribution subtracts 31/(320t^4) inside the functional
multiplied by -hbar/(8pi^2), leaving 4601/(1280t^4). The sign and this
exact coefficient are independently recomputed. Arbitrary finite beta_M
has disappeared because R=0, not because it was chosen zero. Dropping
the negative first-derivative term yields a weaker valid lower bound;
it must never be discarded in the opposite direction. Radiation is only
a metric control, not a claimed solution of the new photon SEE.

## 5. Conditional focusing match and an explicit incompatibility control

For a separately specified actual model obeying

    G_ab+Lambda g_ab=-kappa(T_Maxwell_ab+T_other_ab),

the trace gives R=kappa T_total+4Lambda, hence
R_UU=Lambda-kappa(E_Maxwell+E_other). If E_other>=ell_other, averaging
the field bound gives Q2=kappa B2 and raw constant
kappa B0+Lambda-kappa ell_other. Replacing it by
kappa B0+max(Lambda,0)+kappa max(-ell_other,0) is upwards and nonnegative.
The four source-sign cases are checked independently. Positive radiation
has E=rho>=0. An independent left-hand +alpha I_ab, if present, moves
to +alpha I_ab/kappa in the effective additional source and contributes
-alpha E_I to R_UU; it cannot be silently dropped or conflated with beta_M.

This is a dictionary, not a solution of the SEE. In flat FLRW it is uniform
on comoving normals to constant-time surfaces. It supplies A.1's geometric
QEI hypothesis only when the actual metric/state/source and uniform caps
hold on ALL required windows in its extended domain [-tau/2,tau]. A.1's
initial pointwise bound R_UU<=rho0<=0 on [0,tau/4], the contraction threshold
and its global causal hypotheses still have to be verified separately.
It proves no all-boosted/null bound and removes no initial SEC premise.

If H_j<=c_j/tau^(j+1), D<=tau, A_D is dimensionless and V scales as
tau^-4. With Lambda and other-source penalties displayed separately,
the quantum parts of Q2/tau^2 and Q0*tau^2 scale as kappa hbar/tau^2.
The arithmetic case c_j=1, |beta_M|=1 gives A=23/12,V=864, so these
quantum ratios are 529 kappa hbar/(1152pi^2 tau^2) and
3 kappa hbar/(10pi^2 tau^2). No background or physical prescription is
chosen by those arithmetic controls.

In fact this naive cap has |K|tau=3|H|tau<=3. With zeta=rho0*tau^2=0,
A.1's tail q obeys q(0)=0,q(1)=1 and has nonnegative QEI and L2 terms.
Cauchy--Schwarz gives integral q'^2>=1, so its gradient contribution
(3/(3/4))integral q'^2>=4. Thus that trial cannot focus under |K|tau<=3,
even if the quantum dimensionless ratios tend to zero. The frozen A.1
optimized threshold is likewise >4.4445507. No exclusion is asserted
for other curvature caps, nonzero negative zeta, trial families, surfaces,
or spacetime extensions. This control prevents a small-constant arithmetic
example from being promoted to cosmological incompleteness.
