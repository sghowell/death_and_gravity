# S6.35: exact own-f transfer on a fixed physical-time window

Use the unchanged action, off-shell physical probe g=eta, clock theta=T
and exact positive isotropic own-f branch of S6.33. The physical-g and
clock equations are not imposed. This is a homogeneous linear tensor
theorem with q=0 and specified incoming hidden data, not a full parent
solution or a response excited without data or forcing.

## Uniform actual canonical potential

Write u=T/tau, v=u^2, c=2+delta, k=b^3/N and Y=sqrt(k)Q. The actual
own-f tensor equation is

    Y_uu + V_delta Y = 0,
    V_delta = 2B1 N/b^2 - (sqrt(k))_uu/sqrt(k).

For every |u|<=1/100 and 0<delta<=1/100, the theorem proves

    |V_delta - 16/(delta+8u^2)| < 345169/8000 < 44.

Exact partial Taylor jets and implicit differentiation through third
order give |A_v|<129 and |pump|<22, where A=2Q_profile N/b^2 and
B1=Q_profile/D. A separate univariate Arb jet engine verifies both bounds
on an exact eight-box cover of the entire parameter domain. The positive
root and fixed-point branch are those of the recursively rebuilt parent.

The crucial pole cancellation uses A(0,c)=8c and

    Z=delta+8v, D=Z-v^2 E(v), 0<E<=20, D>=Z/(1+v)^4.

It does not estimate two divergent terms separately. The final strict
rational margin below 44 is 6831/8000. In physical time the potential
remainder is bounded by 44/tau^2.

The regularized background is analytic at its extended corner, but the
pole-subtracted potential is not jointly continuous there: it tends to
24 along v=0 and to 15 along delta=v^2. Both limits use literal positive
delta. There is no literal delta=0 action and no inferred uniform bound
on derivatives of this combined remainder.

## Reference connection and exact endpoint convention

The reference y_xx+16y/(1+8x^2)=0, x=u/sqrt(delta), becomes

    psi_tt + [rho^2+(3/4)sech^2(t)]psi=0, rho=sqrt(7)/2.

This is not the coupled-relative coefficient-80 operator. The right Jost
solution approaches exp(i rho t) on the right and
A exp(i rho t)+B exp(-i rho t) on the left, with

    A=Gamma(1-i rho)Gamma(-i rho)
      /[Gamma(3/2-i rho)Gamma(-1/2-i rho)],
    B=i csch(pi rho), |A|^2-|B|^2=1, |B|>1/40.

Ordinary Gauss connection formulas, Gamma reflection and explicit branch
conventions prove these expressions. The even/odd power connection gives
an independent normalization check. The time-frequency transfer is

    T_t = [[conj(A), -conj(B)], [-B, A]], det(T_t)=1.

The asymptotic radial |u|^(1/2 +/- i rho) connection has determinant -1
because left and right radial orientations differ. It retains the phase
2rho log(4sqrt(2)/sqrt(delta)); it is not itself a finite-endpoint formula.

For the actual finite interval define

    r=sqrt(u^2+delta/8), t=asinh(sqrt(8)u/sqrt(delta)),
    Y=sqrt(r)psi, Z_state=(psi,psi_t/rho).

The endpoint map from (Q,Q_u) is the product

    [[1/sqrt(r),0],[-u/(2rho r^(3/2)),sqrt(r)/rho]]
      [[sqrt(k),0],[k_u/(2sqrt(k)),sqrt(k)]].

The physical input is (Q,Q_T), so its second column additionally contains
tau. The actual k_u is retained. The wave frame

    V_rho(t) = [[exp(i rho t),exp(-i rho t)],
                [i exp(i rho t),-i exp(-i rho t)]]/sqrt(2)

is unitary in Z_state coordinates. With L=1/100 and
t_star=asinh(sqrt(8)L/sqrt(delta)), define exactly

    S_delta=V_rho(t_star)^* F_delta(t_star,-t_star) V_rho(-t_star).

Here F_delta is the actual homogeneous Z_state evolution. A fictitious
free exterior merely supplies these coefficient coordinates. No parent
metric, source, vacuum or quantum state is asserted outside the interval.

## Surviving actual mixing at finite parameter

For every 0<delta<=10^-6 on the fixed physical interval
|T|<=tau/100, the combined theorem proves

    ||S_delta-T_t||_2 < 4/399 < 1/80,
    |B_eff| > 239/15960 > 1/80, B_eff=-[S_delta]21.

The canonical perturbation is E=r^2[V_delta-16/(delta+8u^2)]. Its central
L1 integral is less than 121/25000, and both missing reference tails
together contribute at most 3/3200. Dividing their sum by rho gives
J<2311/500000<1/200. The exact logarithmic norm of the balanced generator
is |W|/(2rho), giving reference evolution norm B0<2. Duhamel's identity
then bounds the coefficient error by 2B0(exp(J/2)-1)<4/399.

The wrapper discharges the connection proof's remainder premise using
the actual coefficient theorem and checks the smaller delta subdomain.
There is no unproved coefficient premise in the combined conclusion.
The error bound need not vanish as delta tends to zero; no complete
outer-dressing limit is claimed.

## Scope boundary

The duration tau/50 is fixed in physical time. This is nevertheless not
a fixed physical low-frequency EFT band, source-matching theorem or
cutoff. Zero hidden data with q=0 still give the zero solution. The
canonical source for nonzero prescribed q is retained algebraically as
r^(3/2)(2B1 b/sqrt(k))q, but a source excitation bound is separate work.

The finite mixing is a specified homogeneous transfer matrix element.
It does not settle the complete covariant CD/matter dictionary, nonlinear
or scalar sectors, controlled all-order reduction, physical vacuum,
positivity/V/G applicability, UV completion or original P8. No retarded
kernel is silently inserted into a single-copy variational action.

The level is CERTIFIED by written continuous proofs, exact arithmetic,
source-audited connection formulas and independent calculations, not
FORMALIZED by a proof assistant. Original P8 remains OPEN.
