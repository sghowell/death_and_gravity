# Full rolling constraints and the decisive matter-cone obstruction

## 1. Actual source and quadratic action

Apply the T=V+U map to the complete 64-component rolling variation,
including every coefficient jet. The stationary background is zero,
and Tstar_0=-d*n, Tstar_i=0, with d=15H/(4h). No lapse-time derivative,
metric shear or vector shift is omitted. The centered mass update
preserves this source and the background; its retained mass matrix on
the clock is diag(1,-1,-1,-1).

For a longitudinal T_i=partial_i sigma, retain its temporal component t,
the original lapse n, scalar shift b and free matter s. The exact added
quadratic density divided by a³ is

    (t+d*n)²/2-q*sigma²/2+zeta*q*(sigma_dot-t)²/2.

The complete original CD/M1 scalar action is present. Eliminating t
gives C*(sigma_dot+d*n)²-q*sigma²/2, where
C=zeta*q/[2(1+zeta*q)]. The shift still fixes
n=(v_dot+ell*s/2)/Theta for finite u!=0. The full auxiliary determinant
is -4q²Theta²(1+zeta*q). Joint solves of all three auxiliary equations
independently reproduce the three-dimensional physical velocity form

    (s_dot+w*v_dot/Theta)²/2
       +J*v_dot²/Theta²+C*(sigma_dot+d*v_dot/Theta)².

Its pivots are 1/2, J/Theta² and C, all positive for zeta>0 and q>0.
Unlike the unchanged-mass trace examples, this candidate has no negative
clock pivot on the actual rolling trajectory.

## 2. Positive principal spatial form is a weaker test than causality

Let f=d/Theta=15/[4(1+u²)^3-1], which is smooth at the center with
f(0)=5 and f_dot(0)=0. In the scalar chart use rho=sigma+f*v and

    hvec=(f_dot,-f*ell/2), x=(v,s,rho), m=(-f,0,1).

The time dependence is retained: the entire extra action becomes
C*(rho_dot-hvec.(v,s))²-q*(rho-f*v)²/2. The old spatial time boundary
has already been included in G0. No derivative of f is dropped and
no new coefficient-dependent time boundary is silently assumed zero.

The new kinetic matrix is diag(K0,C). Its position matrix is

    G = diag(G0-C*hvec*hvec^T/q,0)+m*m^T/2.

All velocity-position terms remain in the action. The Legendre transform
is the exact square (P-linear)^T K^-1(P-linear)/4+q*x^T G*x.
The positive rank-one mass term supplies the third spatial direction.
An exact even-polynomial proof gives

    (hvec^T G0^-1 hvec)/2 < 1000

for every finite real time after the removable factors are cancelled.
The denominator and the degree-50 numerator of the bound are strictly
positive by nonnegative coefficients and positive constant. Since
C/q<=zeta/2, choosing 0<zeta<=1/2000 retains at least half of the
original two-dimensional positive gradient form after the mass square
is minimized over rho. Thus K and G are positive on the punctured
rolling chart for every q>0 in this conservative coupling range.

This is not the original requirement K-G>=0. Positive kinetic and
gradient energy alone does not bound the matter-relative wave cone.

## 3. Exact high-frequency characteristic obstruction

An important exact control is G0=K0 for this original CD/M1 witness:
both of its scalar modes saturate the matter cone. In the limit
q->infinity at a fixed finite u!=0 and nonzero zeta, C->1/2. Lower
position/velocity terms are subleading in the characteristic symbol.
After all scalar constraints and the displayed field transformation,

    det(G_infinity-c²*K_infinity)/det(K_infinity)
       =(1-c²)*[(1-c²)²-kappa*c²],
    kappa=d²/(2J)>0.

The three generalized eigenvalues are exactly

    1,
    c_plus²=1+kappa/2+sqrt(kappa²+4kappa)/2 > 1,
    c_minus²=1+kappa/2-sqrt(kappa²+4kappa)/2 in (0,1).

The latter two have product one. Independently, the (v,rho) principal
minor of G_infinity-K_infinity has determinant -f²/4<0. Both checks use
physical scalars after the full lapse/shift/temporal-vector solve.
The mechanism persists in the same quadratic action with any positive
isotropic retained mass response gamma, replacing kappa by d²/(2gamma J)
and the minor by -f²/(4gamma²). This lemma does not classify every
nonlinear mass-retuned parent.

The result is a high-frequency characteristic limit, not the misleading
finite-q phase speed of a massive particle. Zeta cancels from this
principal limit. Zeta=0 is a different rank limit with no new vector;
one cannot interchange it with q->infinity. Kappa tends to zero at
the center, but the strict violation at every finite u!=0 is enough
to fail the original all-time matter-cone acceptance condition.

Therefore this literal isotropic mass-retuned action is not an accepted
P8 witness, despite exact zero-curl matching and positive principal
energy. No assertion is made that the superluminal regime lies below
a justified EFT cutoff, or that this excludes all other kinetic
operators, source couplings or parents. It does not close original P8.

## 4. First-order crossing and transverse controls

Let pi be the volume-normalized momentum conjugate to sigma. Perform
the Legendre transformation before solving the lapse. The vector
temporal equation gives t=pi-d*n. With the original R0 and C0,

    n=-(R0+d*pi/2)/J,
    H=C0+(R0+d*pi/2)²/J+(1+1/(zeta*q))*pi²/2+q*sigma²/2.

This has no inverse Theta. The actual canonical momenta are
P_b=2a³qv, P_s=a³pm and P_sigma=a³pi; the moving old boundary
-Hubble*b*P_b is retained. At u=0, d=0, so the center momentum Hessian
splits into the old two-dimensional block and the positive vector
block. The regular center velocity chart is positive for q>6; the
first-order coefficients are smooth for every q>0. The new derivative
d_dot(0)=15 is explicitly retained, not declared equal to the old jet.

Both transverse vector polarizations decouple at quadratic order
because Tstar has no transverse source and its background is zero.
For one coordinate polarization their full FLRW action is
a*[zeta*T_dot²-(zeta*q+1)*T²]/2. Its Hamiltonian is
P²/(2a*zeta)+a*(zeta*q+1)*T²/2, positive for zeta>0. The original
metric tensor sector is unchanged. These checks do not undo the
longitudinal/matter-cone obstruction.

## 5. Units and remaining obligations

After extracting M² tau², zeta=zeta_physical/(M² tau²). The retained
isotropic mass parameter is m_physical²=M²/zeta_physical. The nonunit
control M²=3, tau=2, zeta_physical=1/1000 gives zeta=1/12000,
m_physical²=3000 and m_normalized²=12000. These are isolated Proca
parameters, not a certified coupled gap, matching cutoff or UV scale.

A different candidate must address the cone obstruction explicitly.
Neither modifying the mass while retaining this temporal mixing, nor
making the curl arbitrarily small, supplies a literal causal-parent
verdict here. Full nonlinear secondary constraints, causal matching
with quantitative remainders, vacuum/finite-gravity V/G/B estimates
and original P8 closure remain open research, not user-choice blockers.
