# Phase-resolved all-sampler proof

This proof uses the actual positive prepared reference state of A.8. The
auxiliary flat kernel introduced below is only a norm-comparison device;
its modes are not asserted to solve the perturbed field equation.

## 1. Dimensionless variables and the difference inequality

Keep A.8's metric, free real massless minimally coupled field, Hadamard
preparation and named renormalization scheme. Put

    x=eta/eta_star, s=t/(A*eta_star²)=y²/2,
    a_tilde=a/(A*eta_star), ds=a_tilde dx,
    k=k_comoving*eta_star, alpha=alpha_conformal*eta_star.

The unchanged past has x=y. On the history ending at y=3 the potential
U here means eta_star² times the physical conformal potential. The
dimensionless duration from preparation at x=1 is at most T=3. A.7 gives
`|U^(j)|<=delta*b_j`, including j=0,1,2,3, and all past jets vanish.

For a real smooth proper-time sampler h with compact support in 2<y<3,
write `psi(s)=h(A*eta_star²*s)` and `F=a_tilde^(-3/2)*psi`. Let w_k be
sqrt(2*k_comoving) times the actual normalized rescaled field mode, with
past value exp(-ikx), so `w_k''+(k²+U)w_k=0`. Also set
`h_c=a_tilde'/a_tilde` and `D=partial_x-h_c`.
Multiplying A.7's modes by the constant phase needed for this past
convention changes neither the state nor any spectral squared norm.

For an arbitrary Hadamard target state omega (not necessarily homogeneous,
quasifree, zero-mean or equal to the reference omega_delta), the difference
of massless minimal effective energy densities is the state difference
of the square of the proper-time field derivative. The positive-type
worldline argument, applied after the exact rescaling, gives

    integral h²(E_omega-E_delta) dt >= -Q_delta[h],
    Q_delta = hbar/(4*pi³*(A*eta_star²)³) * I_delta[F],
    I_delta = integral_0^infinity k dk integral_0^infinity d_alpha
                  |integral exp(-i*alpha*x) F(x) D w_k(x) dx|².

The argument uses positivity of each differentiated target two-point
function, the common commutator, and the smooth symmetric Hadamard state
difference. It is a statement about the pulled-back positive-type kernel,
not a substitution of conformal time into a proper-time geometric theorem.
The coefficient follows from the reference measure hbar*k dk/(4*pi²)
and the half-Fourier diagonal formula 1/pi. The finite bound below also
justifies the required reference positive-frequency integral. Target states
belong to the same field algebra on the prepared spacetime, so the common
commutator is retained. The reference is prepared on its unchanged past.

## 2. Uniform phase-resolved mode derivatives

Write `w_k=exp(-ikx)*b_k`. The exact equation and initial data give

    b_k''-2ik*b_k'+U*b_k=0,
    b_k'= -integral_1^x exp(2ik(x-r))*U(r)*b_k(r) dr,
    b_k-1= -integral_1^x [(exp(2ik(x-r))-1)/(2ik)]*U(r)*b_k(r) dr.

The second kernel is bounded both by x-r and by 1/k. Ordered Volterra
iteration therefore gives

    |b_k| <= cosh(sqrt(delta*b0)*T)
            <= exp(delta*b0*T²/2) < 2

on A.7's closed amplitude interval. Indeed `(2n)!>=2^n*n!` proves the
middle bound termwise, and `delta_bar*b0*T²<=1` gives the last bound.
At delta=0 the actual solution is b_k=1. The strict bound by 2 is merely
a convenient uniform majorant, not a claim of nonzero perturbation.

For `I[q](x)=integral_1^x exp(2ik(x-r))*q(r)dr`,

    (I[q])'=I[q']+exp(2ik(x-1))*q(1).

Every past derivative of U*b_k vanishes. Hence, exactly,
`b_k^(j)=-I[(U*b_k)^(j-1)]` for j=1,2,3 (and higher when needed).
This step prevents spurious powers of k from differentiating the phase.
It is false without the stated past jets; that omission is a negative
control. The derivatives are of finite smooth integrals before any mode
integration is interchanged.

Let `|b_k^(j)|<=delta*n_j`, j=1,2,3. Direct product bounds give

    n1=2*T*b0,
    n2=T*(2*b1+delta_bar*b0*n1),
    n3=T*(2*b2+delta_bar*(2*b1*n1+b0*n2)).

Define, for j=0,1,2,3,

    N_j=2*b_j+delta_bar*sum_{ell=0}^{j-1}
                         binomial(j,ell)*b_ell*n_(j-ell).

Then `||(U*b_k)^(j)||<=delta*N_j`. One integration by parts in I gives
`|I[q]|<=(||q||+T*||q'||)/(2k)` when q(1)=0. It follows that

    |b_k-1|<=delta*l0/k, l0=2*T*b0,
    |b_k^(j)|<=delta*lj/k,
    lj=(N_(j-1)+T*N_j)/2, j=1,2,3.

The l0 bound uses the second exact Volterra identity directly. No expansion
in delta or assumption k²+U>0 is used; the estimates cover every k>0.
In particular the infrared part is not replaced by an adiabatic vacuum.

## 3. Two derivatives suffice for the positive-frequency error

Split the actual differentiated mode into

    D w_k=exp(-ikx)*[(-ik-h_c)+e_k],
    e_k=-ik*(b_k-1)+b_k'-h_c*(b_k-1).

A.7's `|h_c|<=2`, the identity `h_c'=-U-h_c²`, and its derivative
give caps `H0=2`, `H1=4+delta_bar*b0`,
`H2=delta_bar*b1+4*H1`. On k>=1 the preceding inverse-frequency
bounds imply, for j=0,1,2,

    ||e_k^(j)|| <= delta*E_j,
    E_j=l_j+l_(j+1)+sum_{ell=0}^j binomial(j,ell)*H_ell*l_(j-ell).

Let the conformal sampler support length be ell<=T=3. Two sampler
integrations by parts have no boundary terms and yield

    |Fourier(F*e_k)(alpha+k)|
       <= delta*sqrt(ell)*B/(alpha+k)²,
    B=E0*||F''||+2*E1*||F'||+E2*||F||.

All norms here are in dx. The exact ultraviolet moment is

    integral_1^infinity k dk integral_0^infinity (alpha+k)^(-4) d_alpha
        = integral_1^infinity 1/(3k²) dk = 1/3.

Thus `I_error,UV <= delta²*ell*B²/3 <= delta²*B²`. This explicitly
controls both frequencies; a coincident Wick-square bound would not do so.
It needs only F in H2_0, not a third sampler derivative.

For 0<k<=1 use the nonsingular Volterra kernel instead:
`|b_k-1|<=delta*T²*b0=9*delta*b0`,
`|b_k'|<=6*delta*b0`, so `|e_k|<=33*delta*b0=:delta*E_IR`.
The product F*e_k is generally complex. FULL Parseval, not an unjustified
half-frequency equality, gives

    integral_0^infinity |Fourier(F*e_k)(alpha+k)|² d_alpha
        <= 2*pi*delta²*E_IR²*||F||².

Since integral_0^1 k dk=1/2,
`I_error,IR <= pi*delta²*E_IR²*||F||²`. The disjoint IR and UV
bounds include k=1 by harmless endpoint overlap of measure zero.

## 4. Auxiliary flat kernel and exact clock conversion

Let I0 be the same spectral functional with b_k replaced by 1. The
real function h_c is still the ACTUAL h_c. Integrating k from 0 to
u=alpha+k and applying Parseval gives exactly

    I0=(pi/4)*integral [F''²+2*((h_c*F)')²
                                      +(8/3)*F''*(h_c*F)'] dx.

The density is `(F''+(4/3)*(h_c*F)')²+(2/9)*((h_c*F)')²`.
Thus I0 is nonnegative, even though these auxiliary modes are not a
perturbed-metric state. Cauchy--Schwarz and sqrt(2)<3/2 imply

    sqrt(4*I0/pi) <= ||F''||+(3/2)*||(h_c*F)'||.

The spectral triangle inequality applies to the genuine norm in
L2(k dk d_alpha), once each of the two pieces is finite:
`sqrt(I_delta)<=sqrt(I0)+sqrt(I_error)`.

For the proper variable s, set H=a_tilde_s/a_tilde and a dot to mean
d/ds. The exact transformed integrands, including the square-root
Jacobian dx=ds/a_tilde, are

    F/sqrt(a_tilde)=a_tilde^(-2)*psi,
    F'/sqrt(a_tilde)=a_tilde^(-1)*(psi_dot-(3/2)*H*psi),
    F''/sqrt(a_tilde)=psi_ddot-2H*psi_dot
                                  +[(3/4)H²-(3/2)Hdot]*psi,
    (h_c*F)'/sqrt(a_tilde)=H*psi_dot+(Hdot-H²/2)*psi.

On the actual target `a_tilde=(4s²-delta)^(1/4)`, 2<=s<=9/2.
Writing z=delta/(4s²)<=1/32 gives
`H=1/(2s(1-z))`, `Hdot=-(1+z)/(2s²(1-z)²)`. Consequently

    H<=8/31, |Hdot|<=132/961, a_tilde^(-2)<=1/2,
    a_tilde^(-1)<=1.

These hold on the larger domain delta<=1/2. For example
a_tilde^4>=16-1/2>4 proves the inverse-square cap. On every compact
subinterval of [2,9/2], the two Dirichlet Poincare inequalities and pi>3
give `||psi_dot||<=(5/6)||psi_ddot||` and
`||psi||<=(25/36)||psi_ddot||`. The constants use the full envelope
length, so they hold uniformly for arbitrarily small sampler supports.

Write D2=||psi_ddot||. The exact rational majorants are

    ||F||<=f0*D2,       f0=25/72,
    ||F'||<=f1*D2,      f1=205/186,
    ||F''||<=f2*D2,     f2=9271/5766,
    ||(h_c*F)'||<=b1s*D2,
    b1s=(8/31)*(5/6)+(132/961+32/961)*(25/36).

Let `C0root=f2+(3/2)*b1s`, and define the safe rational error root

    Rroot=2*(E0*f2+2*E1*f1+E2*f0)+2*E_IR*f0.

The ultraviolet factor is actually 2/sqrt(pi)<2; replacing it by 2
and using sqrt(A+B)<=sqrt(A)+sqrt(B) gives this Rroot without hiding
a Fourier or support-length factor. The result is

    I_delta <= (pi/4)*(C0root+delta*Rroot)²*||psi_ddot||².

Exact arithmetic gives `C0root<2.109`, `Rroot<64466343`, and
`(C0root+10^-14*Rroot)²<5`. These are derived from the pinned A.7
potential jets, not supplied mode-error assumptions or floating samples.
All polynomials used in the amplitude majorants have nonnegative
coefficients; the endpoint check controls the entire closed interval.

Since `||h_ddot||²_dt=(A*eta_star²)^(-3)*||psi_ddot||²_ds`, all
dimensionful factors cancel exactly as required:

    Q_delta[h] <= 5*hbar/(16*pi²)*||h_ddot||²_dt,
    0<=delta<=10^-14.

## 5. Absolute inequality and scope

A.8 proves E_delta>0 throughout this target in the SAME named raw-H_lambda,
gamma=0 prescription and amplitude regime. Dropping its positive credit
in the difference inequality proves

    integral h(t)² * [<T_uu>_omega-<T>_omega/2] dt
        >= -5*hbar/(16*pi²)*integral |h_ddot(t)|² dt

for every admitted Hadamard target state and every real smooth compact
sampler in 2<y<3. The renormalization-independent difference inequality
alone would not justify dropping an unknown reference term. Changing
finite curvature coefficients can change the absolute bound; no such
change is covered silently. The coefficient 5 is sufficient, not sharp.

Density extends the statement to real H2_0 samplers on any such compact
interval: the scale/clock maps are smooth and nonsingular, the reference
functional is H2-bounded, and each fixed Hadamard state's renormalized
EED is smooth on the compact support. No bound on the target state's
pointwise EED or energy is assumed. Spatial homogeneity makes the same
reference bounds valid along every comoving line, without restricting
the target state to be homogeneous.

This is a field inequality on the specified off-shell prepared geometry.
It does not turn A.8's small SEE residual into a self-consistent solution,
control any nearby metric not satisfying these bounds, supply the entire
geodesic domain required by A.1, remove A.1's initial pointwise curvature
assumption, cover massive/interacting fields, or prove cosmological
incompleteness. Those remaining P8(a) obligations are not closed here.
