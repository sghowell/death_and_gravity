# Vector determinant, local one-loop potential and scope

The classical source centering does not alter the W Hessian, as the
S6.46 completion of the square established. At fixed light fields
rescale W_c=sqrt(zeta_physical) W. In the physical timelike rest frame
the canonical temporal and spatial mass coefficients are a*m0^2 and
b*m0^2, where m0^2=M^2/zeta_physical, a=1/gamma_t and b=1/gamma_s.
The inherited closed-tube bounds imply

    8/9<a<19/18,    35/36<b<19/18.

These positive coefficients permit the usual constant-coefficient
massive-vector Wick prescription, with a positive Euclidean quadratic
form. This is the local term with curvature and derivatives of the
mass tensor set to zero for its calculation. The coefficients still
depend on phi,X. Freezing them does not construct a global solution
with constant phi and simultaneously nonzero scalar gradient, or a
Lorentz-invariant vacuum. It is not a vacuum scattering calculation.

## 1. Retain the temporal constraint and all polarizations

By spatial rotation take p=(omega,k,0,0). The full Euclidean kernel is

    K_E=(omega^2+k^2)I-p p^T+m0^2 diag(a,b,b,b).

The two transverse entries give (omega^2+k^2+b*m0^2)^2. The coupled
temporal/spatial-longitudinal determinant is

    a*b*m0^4+a*m0^2*omega^2+b*m0^2*k^2
      =a*m0^2[omega^2+(b/a)k^2+b*m0^2].

The full four-by-four determinant and an independent non-aligned
spatial-momentum/rank-one determinant lemma check this exactly.
At a=b=1 it is m0^2(p^2+m0^2)^3, not a two-mode Maxwell determinant.
At k=0 there are still three massive polarizations.

An independent Lorentzian constraint check starts from

    L=a*m0^2*t^2/2-b*m0^2*k^2*sigma^2/2
       +k^2(dot(sigma)-t)^2/2.

Eliminating t gives kinetic coefficient
a*m0^2*k^2/(a*m0^2+k^2) and frequency squared
b*m0^2+(b/a)k^2, agreeing with the Euclidean factor. The frozen-light
longitudinal ratio b/a lies in (35/38,19/16). It is not being used
as a causal-cone verdict for the fully coupled metric/clock system.

## 2. Dimensional regularization and the finite vector constant

Set d=4-2epsilon_DR. Continue the spatial rotation group to d-1
dimensions: there are d-2 transverse polarizations and one longitudinal
polarization. The momentum-independent factor a*m0^2 contributes an
integral of a constant. It is scaleless and vanishes in dimensional
regularization. This statement does not delete a field-dependent
contact measure in another regulator; nor does it account for all
affine-complement or light-field determinants.

Let c^2=b/a and m_s^2=b*m0^2. The spatial longitudinal momentum
rescaling k_new=c k contributes c^(-(d-1)). The complete momentum
dependent determinant therefore has weight

    (2-2epsilon_DR)+(a/b)^(3/2)*exp[epsilon_DR log(b/a)]

multiplying the single scalar logarithmic integral at mass m_s.
For one real scalar that integral is

    -mu_DR^(2epsilon_DR)*(m_s^2)^(d/2)*Gamma(-d/2)
       /[2(4pi)^(d/2)].

Use the modified minimal-subtraction scale
mu_DR^2=mu^2*exp(EulerGamma)/(4pi). Dividing by m_s^4/(64pi^2),
the integral is -2 exp[epsilon_DR(EulerGamma-log(m_s^2/mu^2))]
Gamma(epsilon_DR)/[(epsilon_DR-1)(epsilon_DR-2)]. Its exact pole
and finite expansion is

    -1/epsilon_DR+log(m_s^2/mu^2)-3/2+O(epsilon_DR).

The code checks this Gamma-function normalization and pole residue,
not only an assumed scalar finite constant. Multiply by the full
dimensional vector weight before subtracting the pole. With
L=log(m0^2/mu^2), the result is

    C=2b^2+a^(3/2)b^(1/2),
    V_pole=-m0^4*C/(64pi^2 epsilon_DR),
    B=2b^2[L+log(b)-1/2]
       +a^(3/2)b^(1/2)[L+log(a)-3/2],
    V_MSbar=m0^4*B/(64pi^2).

An independent sign control differentiates the Euclidean scalar
1/2 Tr log(p^2+m^2) twice in m^2. With s=p^2 in four dimensions,
this gives -1/(32pi^2) times the integral of s/(s+m^2)^2. Its exact
primitive is log(s+m^2)+m^2/(s+m^2). The ultraviolet logarithm is
therefore negative, agreeing with the defined epsilon_DR pole.
At d=3 the same Gamma integral gives -m^3/(12pi), also checked exactly.
The reference's printed opposite epsilon definition is not silently
identified with epsilon_DR; see the limited [comparison](literature.md).

In particular a=b=1 gives C=3 and B=3L-5/2, or the standard
massive-vector finite constant 5/6 inside its three-mode factor.
Replacing the vector by three four-dimensional scalar determinants
would instead give B=3L-9/2. The missing dimensional transverse
multiplicity changes the finite answer by 2 and is explicitly tested.
Isotropic a=b at several nonunit masses is a further independent control.

Also dB/dL=C, so the vector contribution to mu*dV_MSbar/dmu is
-m0^4*C/(32pi^2) when these background coefficients are held fixed.
Physical scale independence requires the corresponding running and
matching counterterms. The sign of this scheme-dependent finite
potential is not a vacuum or UV positivity test.

## 3. Actual clock dependence is not removed by source alignment

Insert the original mass functions a(p_affine),b(p_affine), where
p_affine^2=(h-1+N^-2)/(4h). At N=1,

    p_affine=1/2, p_affine,N=-1/(2h),
    p_affine,NN=3/(2h)-1/(2h^2).

At the fixed physical scale mu=m0, the exact coefficient jets are

    C=3,    C_N=20/(9h),
    C_NN=-4(3645h-236)/(2187h^2),
    B=-5/2, B_N=-22/(27h),
    B_NN=2(8019h+4784)/(6561h^2).

The chain rule through p_affine and independent direct N differentiation
of the original mass functions agree. The first derivatives do not
vanish. Thus a zero classical source or a zero first source variation
does not make the determinant independent of the clock fluctuations.
These are derivatives of the local coefficient at fixed u, not the
full quantum lapse/constraint Jacobian. The physical volume, curvature,
mass-derivative terms, canonical changes and renormalization conditions
would all enter that latter calculation.

## 4. A continuous bound on this local term

On the original tube, C<=3(19/18)^2<4. Set mu=m0. From
log(z)=integral_1^z dt/t and the positive a,b intervals,

    |log(a)|<=1/8,    |log(b)|<1/16.

Since a^(3/2)b^(1/2) and b^2 are at most (19/18)^2,

    |B|<=[2(9/16)+13/8](19/18)^2=3971/1296<4.

Consequently |V_MSbar|<m0^4/(16pi^2). Define Lp=M*tau and R=m0*tau.
Relative to the reference density scale M^2/tau^2, this is bounded by
R^4/(16pi^2 Lp^2)<R^4/(144 Lp^2). The strict pi>3 bound follows from
an inscribed regular hexagon, not a floating-point approximation.

For Lp=10^12 and R=1000, m0/M=10^-9,
zeta=zeta_physical/(M^2 tau^2)=10^-6 and zeta_physical=10^18.
The bound is 1/(144*10^12)<10^-14. This zeta also lies in the prior
background-frequency and prepared-response domains. It establishes
smallness of this scheme-specific local term relative to the named
reference scale only. The reference is not a positive lower bound
on the actual energy or an observable error at the bounce.

## 5. What has not been calculated

No curvature, mass-derivative, nonlocal/in-in or higher-loop remainder
has been bounded. No full quantum stress tensor, renormalized bounce,
light/heavy spectrum or interacting cutoff has been established.
The derivative expansion has not been assigned a numerical remainder
from this leading coefficient. Finite matching counterterms, the affine
complement/measure and light/gravity loops remain separate.

The old light cones exactly saturate the physical matter cone. Merely
making a correction small in norm cannot prove preservation of a
semidefinite condition with zero margin. Its structure/sign, or a
separately named controlled matching deformation, still needs analysis.
This is the original cone requirement, not a new completion criterion.
No frozen action is modified, and original P8 remains open.
