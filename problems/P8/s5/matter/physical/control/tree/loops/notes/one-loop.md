# Isolated free-matter counterterm: derivation and limits

## 1. Source input and coefficient convention

For one real scalar, use GammaE=(hbar/2) Tr log(DE/mu²),
DE=-nablaE²+xi RE. The general Laplace-type coefficient, with E=-xi R
and zero bundle curvature, gives

    A2(xi) = (Riem²-Ric²)/180 + (xi-1/6)² R²/2
             + (1/30-xi/6) BoxR.

Here A2 multiplies s² after stripping (4 pi s)^(-n/2); the source calls the
integrated coefficient a4. The BoxR representative is taken directly from
equation 4.28, not its spin table. Local Lorentzian counterterms use the
same curvature polynomial; no global Euclidean CD determinant is needed.
[Vassilevich, equations 1.16, 1.18, 3.6, 4.28 and section 2.3](https://arxiv.org/pdf/hep-th/0306138).

The frozen M1 field is massless and minimal: xi=0. Its background
chi_dot=M/[10 tau (1+u²)^6] obeys (a³ chi_dot)'=0. Its quadratic fluctuation
operator on a fixed metric is independent of that background amplitude;
there is no scalar potential or chi self-interaction in this determinant.
Introducing a positive auxiliary infrared regulator in the proper-time
residue test changes neither the UV A2 residue nor the massless M1 action.
It is not a massive derivative expansion or a computation of infrared
physics.

Define in four dimensions

    C² = Riem² - 2 Ric² + R²/3,
    E4 = Riem² - 4 Ric² + R².

Direct algebra then yields the minimal coefficient

    A2 = C²/120 - E4/360 + R²/72 + BoxR/30
       = E4/180 + Ric²/60 + R²/120 + BoxR/30.

Modulo Euler and total derivatives, the bulk representatives are
C²/120+R²/72 and Ric²/60+R²/120. The conformal-coupling control xi=1/6
removes R²/72 but leaves C²/120; it is not a replacement for M1.
R=0 alone also does not remove C². Conversely, on a Ricci-flat metric the
coefficient reduces to E4/180: a useful check of the four-dimensional
on-shell/topological redundancy, not a CD equation of motion.

## 2. Pole, Wick and scale bookkeeping

Take n=4-2 epsilon. The proper-time A2 contribution is

    -hbar/[2 (4 pi)^(2-epsilon)] mu^(2 epsilon)
      Gamma(epsilon) r_IR^(-epsilon) integral sqrt(gE) A2.

Since epsilon Gamma(epsilon) tends to one, the Euclidean pole is
-hbar integral A2/[2 (4 pi)² epsilon], equivalently
+hbar integral A2/[(4 pi)²(n-4)]. The counterterm has the opposite sign.
The scalar is real: forgetting the determinant's one-half doubles it.

P8 uses +---, R=-6(Hdot+2H²), and the action -M²R/2. Our local continuation
is dt=-i dtE, gL|continued=-gE, GammaE=-i GammaL|continued. Consequently
R_L|continued=-R_E; R², C² and BoxR acquire no remaining metric-sign
factor. For their action coefficients cE=-cL.

| Quantity per integrated A2 | Euclidean | P8 Lorentzian |
| --- | --- | --- |
| Effective-action pole, coefficient of 1/epsilon | -hbar/[2(4 pi)²] | +hbar/[2(4 pi)²] |
| Counterterm, coefficient of 1/epsilon | +hbar/[2(4 pi)²] | -hbar/[2(4 pi)²] |
| Effective-action derivative in ln mu | -hbar/(4 pi)² | +hbar/(4 pi)² |
| Compensating local beta, per A2 coefficient | +hbar/(4 pi)² | -hbar/(4 pi)² |

The scale row follows from GammaE_ren=-(hbar/2) zeta'(0)
-(hbar/2) ln(mu²) zeta(0), with its local UV coefficient identified as
integral A2/(4 pi)². Finite normalization data are separate.
[Vassilevich, equation 2.32](https://arxiv.org/pdf/hep-th/0306138).

Only **constant** scales are compared. Writing L=|ln(mu/mu0)| gives the
isolated magnitudes

    |Delta cR| = hbar L/[72(4 pi)²],
    |Delta cC| = hbar L/[120(4 pi)²].

These are not measurements or bounds on cR(mu0), cC(mu0). The other pieces
of the effective action compensate local running, so Delta c is not an
observable by itself. Choosing mu=mu(t) in a local action would introduce
derivatives of its coefficients and invalidate the constant-coefficient
variation formulas below. No such substitution is made.

The strictly four-dimensional Euler integral has zero compact-support bulk
variation, but a dimensionally continued Euler pole cannot be discarded
before renormalization: its evanescent part can leave finite/anomaly terms.
BoxR is a total derivative in this bulk calculation, not an absent boundary
obligation. The finite-time windows of S5.8 are not physical boundaries
endowed here with quantum boundary conditions. This audit computes no
boundary/initial-surface coefficient. Anomaly and nonlocal terms require
additional information beyond this local-source bound.
[Vassilevich, sections 7.1 and 8.2](https://arxiv.org/pdf/hep-th/0306138).

## 3. Off-shell finite-basis nonclosure, with the frame kept fixed

There is a simple local witness entirely inside the clock tube. Set
phi=t and N=1, hence the clock kinetic invariant is one, but let the flat
FLRW scale factor a(t)>0 be arbitrary. Modulo Euler and BoxR, the isolated
minimal scalar's bulk density becomes

    a³ R²/72 = (a/2) (a_ddot + a_dot²/a)².

Its Hessian with respect to a_ddot is a>0. Its metric Euler-Lagrange
expression is

    a a'''' + 2 a' a''' + (3/2)(a'')²
      - 6(a')² a''/a + (3/2)(a')⁴/a².

In particular the fourth-derivative principal coefficient is a. On this
same restriction, F2(phi,X)R is at most linear in a_ddot, while the terms
quadratic in phi_munu, and canonical matter, depend on at most a_dot.
Arbitrary smooth coefficient functions become functions of t; they do not
alter the zero acceleration Hessian of the strict class. A total derivative
has identically zero Euler-Lagrange expression and cannot remove the
nonzero fourth-order symbol. Therefore the isolated determinant is not
absorbed by any retuning inside the strict finite class in these variables.
No extrapolation to phi=constant or X=0 is used.

This statement is deliberately off shell and frame-specific. It is not
equivalence modulo arbitrary equations of motion or derivative field
redefinitions. The distinction matters: pure Einstein gravity has familiar
on-shell curvature counterterm redundancies, whereas adding scalar matter
changes that conclusion. We borrow no numerical coefficients from that
different theory. ['t Hooft and Veltman, sections 6 and 8](https://www.numdam.org/article/AIHPA_1974__20_1_69_0.pdf).

For an explicit algebraic illustration, consider Einstein plus free chi
only, not the CD DHOST equations. A perturbative redefinition

    delta g^ab = [alpha R^ab + beta R g^ab]/M²

changes -M²R/2 by
-alpha Ric²/2+(alpha/4+beta/2)R². Choosing alpha=1/30, beta=-1/30
trades away Ric²/60+R²/120 from that part of the action. However
delta S_chi=(1/2) integral sqrt(-g) T_chi_ab delta g^ab induces

    R_ab partial^a chi partial^b chi/(60 M²)
      + R (partial chi)²/(120 M²).

Using the illustrative Einstein-scalar equation R_ab=partial_a chi
partial_b chi/M² trades this for (partial chi)^4/(40 M^4). The matter is
no longer free and minimal in the newly named metric. CD also has clock
terms whose variations must be retained. Calling the new variable the
physical metric while dropping these interactions changes the prescribed
theory. A legitimate reparametrization retains the original physical
metric as an observable and carries all induced operators. No complete CD
EOM-reduced operator classification is asserted here.

## 4. Exact CD invariants and metric variation

Put u=t/tau, d=1+u², ell=tau sqrt(d), x=u/sqrt(d), v=x² in [0,1].
Then H=4x/ell and Hdot=(4-8x²)/ell². The compact variable x is not the
clock kinetic invariant. Direct contractions give

    ell² R       = -24(1+6v),
    ell⁴ Ric²    = 192(1+8v+28v²),
    ell⁴ Riem²   = 192(1+4v+20v²),
    ell⁴ E4      = 1536 v(1+2v),          C²=0,
    ell⁴ BoxR    = -240-192v+3456v².

For a quantity ell^(-w) f(x), its cosmic-time derivative is

    d_t[ell^(-w) f(x)]
      = ell^(-w-1)[(1-x²) d_x f - w x f].

The -w x drift is essential; tests independently compare ordinary
cosmic-time derivatives. All these curvature coefficients have their
displayed physical length dimensions. For example |R|<=168/ell²,
Ric²<=7104/ell⁴ and Riem²<=4800/ell⁴. Exact extrema, not sampled times,
establish the bounds.

Vary the covariant action **before** imposing FLRW. Define

    I_R²_ab = (1/sqrt(-g)) delta integral sqrt(-g) R² / delta g^ab
             = 2R R_ab - g_ab R²/2
               + 2(g_ab Box - nabla_a nabla_b)R.

Its lower orthonormal components on CD are

    ell⁴ I_R²_00 = 288(1-16v+36v²),
    ell⁴ I_R²_ii = 576(1-2v-6v²),        i=1,2,3.

The density lies in [-224,6048] (minimum v=2/9); the pressure lies in
[-4032,576]. Trace and conservation give independent algebraic identities
I00-3 Ip=6 BoxR and I00_dot+3H(I00+Ip)=0. More importantly, a separate test
varies the lapse-retaining N a³ R² action with respect to both N and a
before setting N=1. It reproduces both components without using
conservation to define pressure. Dropping Ndot before variation produces
an exact density error 576/tau⁴ at the bounce. The off-shell scale-factor
variation above also matches a² I_R²_ii/12 on CD.

The full C² first variation vanishes when the Weyl tensor vanishes; Euler
and BoxR have zero compact-support bulk variations in four dimensions.
Thus I_A2=I_R²/72 on CD, with component absolute suprema 84/ell⁴ and
56/ell⁴. This is not obtained by declaring all curvature terms absent
after inspecting their background values. Indeed the chosen A2 local
density is (32/3) v(8+37v)/ell⁴ and vanishes at the bounce, while
I_A2_00=4/tau⁴ there. Its pointwise density is not its variational source.

## 5. What the scale estimate actually bounds

For a finite local addition cR integral R², the Einstein-form source
contains 2 cR I_R²_ab (the sign depends on which side of the equation it
is placed). Therefore, componentwise,

    |2 cR I_R²_ab| / (M²/ell²) <= 12096 |cR|/(M ell)².

Use the positive reference scale M²/ell²: the actual G00 vanishes at the
bounce and is not an admissible denominator. For the isolated logarithmic
increment with |log(mu/mu0)|<=L the bound becomes

    168 hbar L/[(4 pi)²(M ell)²]
      <= (7/6) hbar L/(M tau)².

The rational enclosure uses pi>3 and ell>=tau. At the sufficient S5.8
choice M tau=10^324 it is at most (7/6) hbar L 10^-648. This is a
conditional estimate for the named local source, not a bound on the full
renormalized stress tensor or a proof of small displacement of a bounce
solution. Finite matching coefficients remain arbitrary until supplied by
a matching calculation or assumption. Exact finite rational/algebraic
parameters, or symbols with proven sign/finiteness assumptions, are
accepted by the bound interface; rounded floats and infinities are rejected.

## 6. Tensor sensitivity and remaining quantum work

Background C=0 does not set the second variation of integral C² to zero.
The independent static TT test uses h11=f(z), h22=-f(z). Direct
linearized Riemann contraction yields R_linear=0, Ric²=(f'')²/2,
Riem²=2(f'')², E4=0 and C²=(f'')² at quadratic order. Dropping C²
therefore misses a real quadratic coefficient. A separate null-plane-wave
control has vanishing scalar contractions despite nonzero Riemann
components, so the converse is not inferred either.

For a unit-norm flat TT polarization, modulo Euler/boundary terms,

    (integral C²)^(2) = (1/2) integral (Box h_ab)².

With h=2Y/M at a=1 this gives 2 cC (Box Y)²/M². This is a
highest-derivative coefficient identity, not a division by the Einstein
kernel on shell and not a curved finite-band norm estimate. The added
R² term has f_R=-M²/2+2cR R, hence its two-derivative TT coefficient
shift is delta M_T²=-4cR R, bounded in magnitude by
672 |cR|/(M ell)² relative to M². Lower-derivative, Weyl, scalar/matter
and constraint effects are not thereby controlled.

The S5.8 spatial band alone supplies no bound on every time derivative
in the new operators. One still needs matched finite coefficients,
massless nonlocal kernels and a specified Lorentzian state/observable,
the anomaly/finite terms, the other loop sectors, a corrected coupled
constraint reduction, and quantitative perturbative/backreaction errors.
No cancellations between uncomputed sectors are ruled out. No global
vacuum or infrared prescription is furnished by a local heat coefficient.

Finally, treating a truncated curvature-squared action as an exact
higher-derivative theory and identifying its extra resummed poles is not
a proof of physical ghosts in the EFT. Perturbative higher-derivative
equations require a controlled physical-branch/order-reduction treatment.
That treatment is not performed for CD here. [Parker and Simon, sections I–II](https://arxiv.org/pdf/gr-qc/9211002).

S5.9 supplies the first mandatory isolated-loop operator audit and a
conditional local-source estimate. It closes neither total loop control,
radiative cone protection, nonlinear stability nor S6 UV matching.
