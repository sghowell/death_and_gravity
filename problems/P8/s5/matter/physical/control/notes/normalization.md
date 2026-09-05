# Exact coupled normalization and a local free-energy theorem

All statements concern the fixed CD/M1 witness and S5.6.CD Hamiltonian.
M=tau=1 may be used to check the algebra, but physical time, the factor
ell=tau*sqrt(1+u^2), and the a^3 measure are restored in the theorem.
The matter metric is unchanged. These calculations are not an additional
principal-speed assumption or a new covariant witness.

## 1. Generic canonical construction

Write a real mode pair with density

    H2 = p^T*A*p/2 + p^T*B*Q + Q^T*C*Q/2,
    action = integral a^3*(p^T*Qdot-H2) dt,

where A and C are symmetric. The complex Fourier version follows by
polarization, pairing opposite momenta; all energies below are Hermitian.
Choose a positive lower-triangular T with T^T*T=alpha=A^-1. Set

    Y=a^(3/2)*T*Q, pi=a^(3/2)*T^-T*p,
    F=3*H*I/2+Tdot*T^-1,
    Bhat=T*B*T^-1+F, Chat=T^-T*C*T^-1.

The symplectic one-form contributes **+pi^T*F*Y** to the transformed
Hamiltonian. Consequently its quadratic terms are pi^T*pi/2 plus
pi^T*Bhat*Y plus Y^T*Chat*Y/2. Put

    S=sym(Bhat), Omega=-skew(Bhat), P=pi+S*Y,
    W=Chat-S^2+[S,Omega]-Sdot.

Since S is symmetric, shifting pi=P-SY is canonical. Its time boundary
contributes -Y^T*Sdot*Y/2, not zero. The resulting Hamiltonian and
Lagrangian are exactly

    Hnorm=P^T*P/2-P^T*Omega*Y+Y^T*W*Y/2,
    Lnorm=|Ydot+Omega*Y|^2/2-Y^T*W*Y/2.

In particular the original momenta required by the interacting theory are

    Q=a^-3/2*T^-1*Y,
    p=a^-3/2*T^T*(P-S*Y).

The time-dependent generating terms here are quadratic. Higher Hamiltonian
orders undergo this linear canonical substitution; replacing p with a
free velocity instead is not its interacting counterpart.

The free equations and positive-energy candidate are

    Ydot=P-Omega*Y, Pdot=-W*Y-Omega*P,
    E=(P^dagger*P+Y^dagger*W*Y)/2,
    Edot=Y^dagger*(Wdot+[Omega,W])*Y/2.

The commutator has the **displayed sign** and is symmetric. There is no
norm-growth term proportional to |P|*|Omega P|: antisymmetry cancels it
exactly. Equivalently solve Odot=-Omega*O with orthogonal initial O and
write Y=O*Z, P=O*Zdot. The connection disappears and the potential becomes
O^T*W*O, with derivative O^T*(Wdot+[Omega,W])*O. This is an exact change
of frame, not adiabatic diagonalization or a choice of quantum vacuum.

`independent.generic_checks` verifies the canonical, boundary, equation,
and energy identities with arbitrary matrix entries. Explicit nonzero
controls delete the volume generator, Sdot, or connection, or reverse
the covariant commutator. Written matrix calculus, not a finite fixture,
justifies application to all admitted times.

## 2. Compact weights and exact arithmetic

Define x=u/sqrt(1+u^2), z=1/q and the compact matter velocity l=ell*chi_dot.
With the prior dimensionless normalization,

    Hbar=4x, Theta_bar=x*(4-(1-x^2)^3),
    Lambda=1-3*(1-x^2)^3/2, l^2=(1-x^2)^11/100,
    wbar=-l*Lambda, Jbar=ell^2*Jphys, J0bar=Jbar+l^2*Lambda^2/2.

Below J, J0 and theta mean these compact coefficients. Pinned S5.5 proves
1/10<J<8. Since |l|<=1/10 and |Lambda|<=1, also 1/10<J0<9.
For a coefficient with physical dimension ell^-n, differentiation at
**fixed comoving momentum** is

    ell^(n+1)*d_t(f/ell^n)=D_n f,
    D_n=(1-x^2)*partial_x+6xz*partial_z-11xl*partial_l-nx.

Thus D_0 q=-6xq, whereas D_2 q=-8xq. The latter is the physical k^2
drift. Also ell*d_t log(a^(3/2))=6x. Neither coefficient is the D-only
background's drift. The differential operator preserves l^2-(1-x^2)^11/100.

Every stored scalar coefficient is represented exactly as

    even(x,z)+l*odd(x,z), with even,odd in QQ(x,z).

The implementation evaluates the expression tree in this quadratic
extension of the rational function field, reducing l^2 by its defining
relation. Multiplication and inversion are exact; there is no series
truncation, interpolation, floating-point root classification, or sampled
high-q extrapolation. Tests compare this arithmetic with an independent
polynomial-remainder route. It avoids a prohibitively large intermediate
multinomial expansion in the gamma chart without changing the mathematics.

Gamma's configuration has a different physical weight from the matter
configuration. With delta=0 in unitary and delta=1 in gamma, define
Lw=diag(ell^-delta,1). Then

    alpha_phys=Lw*alpha_bar*Lw,
    Bphys=ell^-1*Lw^-1*Bbar*Lw,
    Cphys=ell^-2*Lw*Cbar*Lw,
    Tphys=Tbar*Lw.

These are checked against the original cosmic-time Hamiltonian, including
first and second time jets, before point evaluation. Treating the two
gamma scalar coordinates as having the same weight would be incorrect.

## 3. Positive kinetic charts

Use

    Tbar = [[sqrt(d1),0],[chi*sqrt(d2),sqrt(d2)]],
    alpha_bar=[[d1+chi^2*d2,chi*d2],[chi*d2,d2]].

The factors, derived from the full finite-q inverse Hessian, are

| chart | d1 | d2 | chi |
|---|---|---|---|
| unitary | 2J/theta^2 | 1 | -l*Lambda/theta |
| gamma | 2J/R | R/D | l*Lambda^2/R |

Here R=Lambda^2-Jz, D=Lambda^2-J0z. The exact matrices satisfy both
Tbar^T*Tbar=alpha_bar and alpha_bar*Bbar+beta_bar=0. For q>=1000 use
unitary when |x|>=1/9 and gamma when |x|<=1/4. They overlap.

On the unitary exterior, 3<=theta/x<=4 gives

    1/80<=d1<=144, d2=1, |chi|<=3/10.

On the gamma core, (1-x^2)^3>=(15/16)^3, hence Lambda is negative and
1933/8192<=|Lambda|<=1/2; in particular |Lambda|>1/5. Thus

    R>4/125, D>31/1000,
    R>=(4/5)*Lambda^2, D>=(3/4)*Lambda^2,
    4/5<=d1<=500, 1<=d2<=4/3, |chi|<=5/8.

The last chi bound is deliberately loose. D>0 is precisely the full
gamma positive-velocity condition; the chart's bounce pole is q=6.

For either chart, tr(alpha)=d1+d2*(1+chi^2), det(alpha)=d1*d2. The upper
eigenvalue is at most the trace and the lower is at least det/trace.
Substitution of the listed interval bounds, with monotonicity of
d1/(d1+c) for c>0, proves

    (1/100)*I < alpha_bar < 512*I,
    sqrt(d2/d1)<10.

For example det/trace is bounded below by
(1/80)/(1/80+1+9/100)>1/100 on the exterior; in the core the looser
(4/5)/(4/5+(4/3)*(1+25/64)) also exceeds 1/100. The replay checks every
rational implication and chart overlap. These are bounds in local
dimensionless coordinates, not an assertion that dimensionful gamma
alpha has the same uniform eigenvalue bounds in a fixed unscaled basis.

## 4. Exact luminal remainder and retained mixing

For the Legendre-transformed density set beta=-alpha*B and
Gamma=C-B^T*alpha*B. Integration by parts of the symmetric cross term
gives the effective configuration matrix Gamma+beta_s_dot+3H*beta_s.
If s=(delta,0), its compact principal-subtracted remainder is

    Rmass_ij=Gamma_bar_ij+D_(1+s_i+s_j)(beta_s_bar_ij)
             +3Hbar*beta_s_bar_ij-q*alpha_bar_ij.

Exact cancellation removes every high-frequency pole from Rmass. This
is a full finite-q identity, not the replacement of a finite-q matrix
by its principal limit.

Define r=d2/d1, f1=6x+D_0 log(d1)/2-delta*x,
f2=6x+D_0 log(d2)/2, h=D_delta chi, and k=beta_a_bar_12/d2. Then

    Fbar=[[f1,0],[sqrt(r)*h,f2]],
    Kbar=[[0,sqrt(r)*k],[-sqrt(r)*k,0]],
    Omegabar=Kbar-skew(Fbar),
    Omegabar_12=sqrt(r)*(k+h/2).

An independent Lagrangian route gives

    Rhat=Tbar^-T*Rmass*Tbar^-1,
    M=Rhat-Fbar^T*Fbar+Fbar^T*Kbar-Kbar*Fbar
       -D_1(sym(Fbar))-Omegabar^2,
    W=(q*I+M)/ell^2, Omega=Omegabar/ell.

In differentiating off-diagonal entries, D_0 sqrt(r) is retained. The
implementation stores off-diagonals as sqrt(r) times a rational parity
coefficient and differentiates that factor explicitly. It also records
Sbar=sym(Fbar)-Tbar^-T*beta_s_bar*Tbar^-1 for the inverse canonical map.
Sbar is allowed to contain positive powers of q: it is not silently
classified as a bounded mass remainder.

All M and Omegabar entries and the mass covariant derivative

    Cmass=D_2 M+[Omegabar,M]

are bounded rational parity coefficients on the stated compact charts.
The leading q*I commutes with the connection exactly, so

    ell^3*(Wdot+[Omega,W])=-8xq*I+Cmass.

For a tensor polarization with norm E:E, take
Y_T=a^(3/2)*sqrt(E:E)*h_T/2. Direct normalization gives

    M_T=-6-24x^2, D_2 M_T=-36x+96x^3, Omega_T=0.

Both scalar mass eigenvalues tend to -30 and Omegabar tends to zero at
each compact tail; M_T also tends to -30. These compact limits do not
put a fixed-comoving solution in a high-frequency band at infinite time.

## 5. Exact uniform majorants

For each rational parity coefficient the code forms the exact polynomial
LCM of its even/odd denominators. It divides that denominator only by
the following **already proved nonvanishing factors**:

- J and J0, with absolute lower bounds 1/10;
- in unitary, theta/x and x, with lower bounds 3 and 1/9;
- in gamma, R, D and Lambda, with bounds 4/125, 31/1000 and 1/5.

Any nonconstant unexplained remainder causes failure. Each numerator is
bounded by its weighted coefficient L1 norm using |x|<=1 on the exterior,
|x|<=1/4 on the core, 0<=z<=1/1000 and |l|<=1/10. Dividing by the proved
denominator lower bound gives an exact rational majorant. Enlarging the
numerator domain while keeping the actual denominator domain is a valid
one-sided bound. Zeros of a factor outside the claimed chart are not
discarded by cancellation or numerical tolerance.

For symmetric two-by-two matrices, the operator norm is bounded by the
maximum absolute row sum. Thus if b11,b22,b12 bound the two diagonal
entries and the rational off-diagonal factor, use

    C0=max(b11,b22)+10*b12,

and similarly C1 for Cmass. The antisymmetric connection norm is bounded
by 10 times its factor bound. The exact fractions are in the certificate.
For orientation only (not proof inputs), the dominant exterior bounds are
C0 about 1.724*10^13 and C1 about 5.299*10^19; the core bounds are much
smaller. Tensor bounds are C0=30 and C1=132.

Choose the first decimal power at least 1000, 2*C0 and C1 on both scalar
charts and tensors. Exact comparisons give q_star=10^20. Therefore for
q>=q_star,

    q/(2ell^2)*I <= W <= 3q/(2ell^2)*I,
    norm(Wdot+[Omega,W]) <= (8q+C1)/ell^3 <= 9q/ell^3,
    |Edot| <= (9q/(2ell^3))*|Y|^2 <= 18E/ell.

This conservative sufficient result does not rely on the size of an
instantaneous eigenvector derivative or an assumed WKB solution. It is
not an estimate of cubic or quartic coupling strength.

## 6. One-chart physical windows and energy propagation

For any finite center t0, let I0 be |t-t0|<=ell0/100. Since
ell_dot=x and |x|<=1, throughout I0,

    99/100 <= ell/ell0 <= 101/100,
    |x-x0| <= 1/99,
    |log(a/a0)| <= 4/99,
    95/99 <= a/a0 <= 99/95.

The second and third inequalities follow by integrating |xdot|<=1/ell
and |H|<=4/ell. The last follows from e^-v>=1-v and e^v<=1/(1-v)
for 0<=v<1. Choosing gamma when |x0|<=9/50, unitary otherwise, keeps
that single chart inside its covering domain for the **whole** interval.

Since q/q0=(ell/ell0)^2*(a0/a)^2,

    (19/20)^2 <= q/q0 <= ((101/100)*(99/95))^2 < 2.

Thus q0>=2*q_star suffices throughout the window, with slack. For any
s,t in I0 the full-window integral is at most 2/99, so integration of
the energy-rate inequality yields

    |log(E(t)/E(s))| <= 36/99,
    exp(36/99) <= 1/(1-36/99)=11/7 < 2.

Positive W implies E>0 for a nonzero solution. Existence and uniqueness
for this finite-dimensional smooth linear system then give
7/11<=E(t)/E(s)<=11/7. The zero solution satisfies the undivided energy
inequality. Both real and complex solutions obey the same estimate.

The normalization is nonsingular on each selected window. No crossing
of x=0 occurs in a selected unitary window, and gamma windows cross the
bounce without approaching their positive-velocity boundary. This is not
a chart-stitched global energy estimate, infrared theorem, or statement
about one comoving mode across the entire cosmological history.

## 7. Independent checks and acceptance boundary

`independent.py` never imports the compact model or oscillator. It takes
the pinned Hamiltonian's original u=t/tau functions, holds k_com fixed,
and differentiates the full kinetic factor and Hamiltonian before point
evaluation. Its generic canonical proof and known-answer bounce formulas
are separate from the main Lagrangian formula in section 4.

At u=0 the positive gamma domain is q>6 and the two potential diagonals are

    W11=(239800q^2-5994805q-2875202)/(1199*(200q-1199)),
    W22=(200q^3-3603q^2+21620q-43338)/((q-6)*(200q-1199)),
    W12/sqrt(r)=-80q*(199q+5103)/(200q-1199)^2,
    det(W)=(1199q^3-37198q^2+166106q+86676)/(1199*(q-6)),
    Omega=0.

These formulas match the compact route for symbolic q. At q=8 the
kinetic factors are positive but W11 and det(W) are negative. This is
outside the sufficient band and is only a negative control for the
inference 'positive kinetic implies positive oscillator potential.'
An instantaneous indefinite W in time-dependent canonical coordinates
does **not** establish a physical dynamical instability.

The independent nonzero-time fixtures use u=+-9/40 (gamma) and +-3/4
(unitary), as well as the bounce, with k_com^2=10000. All full W, Omega
and S entries agree exactly after restoring ell weights. In particular
the unitary positive-Cholesky connection has physical entry

    Omega12=-sign(Theta)*sqrt(2)*chi_dot*(Theta-H*Lambda)/sqrt(Jphys),

which is negative on both exterior pieces for this witness, not zero.
The fixture momenta are identity checks; they do not satisfy or establish
the much more conservative all-time q_star bound.

All prior sources and certificates remain immutable. The certificate
records parity-rational coefficient digests and denominator/term metadata
instead of reproducing enormous expanded polynomials; the exact source
recomputes them and all bounding fractions. Its read-only replay also
checks the prior S5.6 source chain. The full P8 suite separately replays
the prior gates. The generic calculus and finite-window inequalities are
written proofs; no machine formalization or independent external review
is claimed.

The next gate is a **new** nonexceptional interacting M1 kernel bound
using the canonical map in section 1. The present result supplies its
free propagator estimate but does not establish that gate, loop control,
technical naturalness, vacuum/tube matching, or full P8 completion.
