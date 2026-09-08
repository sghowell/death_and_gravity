# Dimensional local matching: derivation and explicit boundary

The original P8 and
S6.42 action are unchanged. S6.50/51 give an exact-mode comparison
and a finite prescribed subtraction; this work studies its local
matching. The exact local statements below do not supply full finite
matching of the actual clock-mass derivative terms.

## Radial and geometric calculations

Let D=3-2epsilon_DR and omega^2=m^2+k_com^2/a_s^2. The canonical
normalizations are gT^2=a_s^(D-2) and
gL^2=a_s^D*m^2*q/(q+m^2). Their logarithmic rates are
[(D-2)/2]H and [(D-2)/2+z]H. Direct D-dimensional lapse and
scale-factor variation reproduces the energy and pressure weights.
There are D-1 transverse polarizations and one longitudinal mode.
Dimension here is an analytic regulator, not a claimed state on a
noninteger-dimensional Hilbert space.

For an order 2n coefficient sum_j c_j(D) z^j in the common
omega/(4a_s^D) normalization, the radial integral follows from

    integral d^D p/(2pi)^D omega^(1-2n)*z^j
      = m^(D+1-2n)/(4pi)^(D/2)
        *Gamma(j+D/2)/Gamma(D/2)
        *Gamma(n-1/2-D/2)/Gamma(n+j-1/2).

Define G_n(D)=2sqrt(pi) sum_j c_j(D)(D/2)_j/Gamma(n+j-1/2).
Relative to m^(4-2n)/(64pi^2), with ell=log(m^2/mu^2), the
MSbar-normalized radial Laurent terms are

    P_n=(-1)^(2-n) G_n(3)/(2-n)!,
    E_n=-2*(-1)^(2-n) G_n'(3)/(2-n)!,
    radial_finite_n=P_n[Harmonic(2-n)-ell]+E_n.

This component-wise finite part is not automatically the finite
covariantly renormalized stress component. Counterterm variation
must also precede the dimension limit. The code deliberately calls
it radial_finite_part rather than silently labeling it the matched
energy. At order zero it agrees with the frozen S6.47 finite lapse
term including a_N,b_N, and its pole is -(3+3a_N/2+9b_N/2).

The ordinary-Proca control sets a_N=b_N=0. Its energy poles are

    P0=-3,
    P2=-6H^2,
    P4=H'^2-6H^2H'-2HH''.

They agree with direct lapse variations of the frozen curvature
pole action. Use the four-dimensional pole coefficients, but vary
the invariants in D+1 dimensions before expanding. The boxR term
is a divergence under compactly supported local variations; this
does not set its noncompact infinite-time flux to zero.

For a local density N*a_s^D*F, the lapse energy is minus its
Euler derivative in N divided by a_s^D. The implementation retains
N' in R, Ricci^2 and Riemann^2. Pressure of this covariant diagnostic
can also be recovered by its D-dimensional Ward identity. The
radial finite components have the nonzero four-dimensional defects

    order 2: 8HH',
    order 4: (4H/3)[3HH''+6H'^2+H'''].

The epsilon coefficient of the dimensional pole-counterterm variation
removes these defects. The ordinary-Proca matched local energy terms,
with the same respective m^4,m^2,m^0 prefactors, are

    3ell-5/2,
    H^2(6ell-10),
    (ell+2)(6H^2H'+2HH''-H'^2).

An independent local heat-action calculation gives the same answers.
The universal Laplace-type input is the one already used in S6.49,
from [Vassilevich](https://arxiv.org/abs/hep-th/0306138). The dimension
dependence below follows by retaining the identity-bundle trace while
using the explicit curvature traces already checked in that checkpoint.
For the Hodge one-form-minus-scalar ratio in D+1 dimensions,

    a0=D, a2=(D/6-1)R,
    a4=(D-12)R^2/72+(90-D)Ricci^2/180
        +(D-15)Riemann^2/180+(D-5)boxR/30.

The rank derivative of a4 is the minimally coupled scalar a4.
Direct Gamma expansion gives the finite Euclidean action coefficients

    m^4(3ell-5/2), m^2 R(ell-5/3),
    2ell*a4_at_four+4*(partial_D a4)_at_four,

relative to 1/(64pi^2). The Lorentzian local density has the opposite
sign. Varying this finite action independently matches the corrected
ordinary-Proca energy and pressure. Universal local heat inputs are
the same ones retained and checked in S6.49; no finite nonlocal heat
expansion is inferred.

## Independently reconstructed actual clock-mass poles

Let E=H'^2-6H^2H'-2HH'' and R=6(H'+2H^2), with
boxR=-R''-3HR'. The actual dimensionally integrated extra lapse poles
admit the independently checked covariant reconstruction

    extra_P2=-(3H^2+4H')a_N-(9H^2+2H')b_N,
    extra_P4=(b_N-a_N)E/2+(a_N+b_N)boxR/12.

The second-order formula follows by contracting the Proca coincident
Green pole R*g/12-5Ricci/6 against the mass direction
deltaM/m^2=b_N*g+(b_N-a_N)n*n. For the fourth-order term use
the Hodge decomposition G_Proca=G_one_form+gradient G_scalar/m^2.
Only the scalar-gradient piece has a fourth-derivative UV pole.
For a minimally coupled scalar, the regulated equation of motion
gives Z_munu=T_munu+g_munu*box<phi^2>/4. Its pole
<phi^2>=[-m^2+R/6]/(16pi^2 epsilon_DR) and the scalar curvature
stress pole imply, in units 1/(64pi^2 epsilon_DR),

    Z_nn=E-boxR/6, trace Z=boxR/3.

Contracting the same mass direction yields the predicted extra_P4.
The implementation compares both expressions with the independently
integrated actual energy coefficients. The scalar curvature stress
used here is itself obtained by varying its scalar heat coefficient,
not assumed equal to the vector coefficient. Their FLRW lapse values
agree after that variation; the trace and normal projections are
checked separately.

The Green-function argument concerns local dimensional poles. The
minimal scalar action gives
T_munu=Z_munu-g_munu[(nabla phi)^2+m^2 phi^2]/2 and the regulated
equation of motion gives box<phi^2>=2<(nabla phi)^2>+2m^2<phi^2>.
Thus Z=T+g box<phi^2>/4. The coincident contact distribution is a
scaleless trace in the stated dimensional prescription, not a finite
cochain determinant set to zero. Finite trace/EOM anomalies are not
inferred to vanish: only these pole identities are used. The Hodge
operator factorization is the frozen one from S6.49; the scalar
gradient term is retained, rather than dropping the constrained
polarization. Flat Gamma residues and the actual lower-order mass
direction independently check its sign.

Twenty exact reductions also compare the dimension-dependent WKB
and subtraction coefficients at D=3 with the full S6.50/51 functions
on the actual physical clock. Direct Gamma limits check the residues
and finite terms used in the radial continuation. Freezing the
polarization count before integration loses 2+4b_N in the finite
zeroth-order energy, in the common normalization. A separate check
recovers the already frozen S6.47 finite potential lapse jets.

## Finite extension warning

The actual a_N,b_N terms are not an ordinary conserved Proca stress
component; they exchange with the clock equation. The ordinary Ward
test above is a diagnostic, not a condition imposed on that component.

An evanescent covariant continuation of a mass counterterm need not
equal the continuation used by the already frozen constant-potential
MSbar prescription. For example, a spatial trace in D dimensions
and three times its spatial average coincide at D=3 but have different
epsilon derivatives. No such finite shift may silently change S6.47.
Any full matching must explicitly preserve its finite potential and
lapse jets, and derive the finite derivative terms in the same stated
prescription. State admissibility, higher loops and V/G/B remain open.
