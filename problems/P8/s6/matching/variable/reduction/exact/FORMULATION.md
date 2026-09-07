# S6.33: an exact own-f branch with a uniform positive lapse

Keep the unchanged S6.20 parent, both positive Einstein coefficients M^2,
the original physical metric g, canonical clock and free canonical chi.
The same S6.30/S6.32 off-shell unit-volume probe has g=eta and the clock
label theta=T=tau*u. The canonical field is still the frozen phi(theta),
not a replacement clock kinetic term. We solve the homogeneous isotropic
own-f equations only. The physical-g/clock equations are not imposed.

The actual quantitative domain is

    M,tau>0; c=2+delta; 0<delta<=1/100; |u|<=1/100.

The claim is CERTIFIED by a written continuous proof, exact symbolic
identities, independent rational and Arb enclosures. It is not FORMALIZED,
a controlled EFT theorem, or original P8 completion.

## Exact branch and uniqueness

Write f=diag(N^2,-b^2,-b^2,-b^2), with b,N positive. The literal own-f
action, after exposing the covariant-Einstein boundary, is

    Lf=-3M^2 b b_T^2/N-2beta1 N-6beta1 b-2beta4 N b^3.

The lapse and scale Euler expressions obey the undivided consistency
identity in the proof. For any regular solution b in C2 and N in C1
through b(0)=2, the lapse equation first enforces b_T(0)=0. The full
consistency equation then gives an even algebraic relation with nonzero
b-Jacobian 96/[c(c-2)]. Analytic implicit-function uniqueness forces
evenness of b and then N; parity is not an initial assumption. No
velocity or second derivative is divided out at the turning point.

On the stated window this germ has a unique positive regular continuation,
constructed by a desingularized fixed-point equation. Let v=u^2 and

    d=1+v, D=c-2/d^4, J=c d^4(1-7v)+12v,
    Q=32(1-v)/(c d^14), R=d^8 J/[8c(1-v)],
    L=R_v/R, H=D(J_v/J-6/d)-D_v,
    b=[(1-v D zeta)/R]^(1/3), F=2(v zeta H-L)/3,
    Tmap=3b F^2/(2Q), zeta=Tmap.

R here is not a curvature scalar, Q is not a tensor amplitude, and D is
not a time derivative. The positive cube root is fixed. These are exactly
the unchanged profiles beta1=(M^2/tau^2)Q/D and beta4=-beta1 R.

The full continuous box v in [0,1/10000], delta in [0,1/100], zeta in
[11,13] satisfies strict self-mapping and contraction:

    11<Tmap<13, |Tmap_zeta|<1/100, 19/10<b,N<21/10.

All required denominators and the root argument stay positive, while F
and the total db/dv stay negative. The lapse is recovered from the actual
implicit derivative,

    zeta_v=Tmap_v/(1-Tmap_zeta),
    N=2(b_v+b_zeta zeta_v)/F.

Banach existence and analytic-IFT patching establish the branch on the
whole box; the undivided identity proves both Euler equations at and away
from the center. An independent Fraction engine and a separate 256-bit
Arb evaluation certify the entire box, without sampling or fitted roots.

The delta=0 face is only a joint analytic extension of these equations
and metric functions. At delta=v=0 the original action has D=0 and Q=16
and is undefined. The certificate does not add a literal c=2 theory.

## Exact center and tensor equation

The actual center values include

    b0=2, N0=c^2/2, b_uu(0)=-2c(c+2),
    N_uu(0)=(27c^4-12c^3-236c^2+592c-512)/4,
    f00(0)=c^4/4,
    f00(0)-(f0+f2)00(0)=(c^2-4)^2/4.

This is a different off-shell metric component from the original rolling
solution's f00=c^2. It is not a covariant CD coefficient match or a bound
on the complete derivative expansion.

For g=diag(1,-exp(q),-exp(-q),-1) and hidden anisotropy Q_tensor, the
literal quadratic action on this exact own-f background is

    LT=M^2 q_T^2/4+Kf Q_tensor,T^2/4-beta1 b(Q_tensor-q)^2/2,
    Kf=M^2 b^3/N.

A separately derived four-coordinate connection/Ricci calculation
reproduces the action and both tensor Euler equations, including the
Einstein boundary. Keeping q prescribed gives

    D_T(Kf D_T Q_tensor)+2beta1 b(Q_tensor-q)=0.

At the center Kf=16M^2/c^2 and the own-f algebraic mass squared is
8c/[tau^2(c-2)]. On each fixed compact inner interval x=u/sqrt(delta),
joint metric analyticity controls the drift and canonical pump and yields

    Q_tensor,xx+16(Q_tensor-q)/(1+8x^2)=0.

If q is instead varied too, its limiting equation contributes 64 and
the coupled relative mode has coefficient 80/(1+8x^2). These two
operators are distinct. ODE solution convergence requires convergent
inner-coordinate initial data and prescribed inputs. This fixed-inner-
interval statement is not a fixed physical-window response bound.

## Unfinished obligations

No tensor initial data or Green prescription is chosen by this background
construction. It supplies neither an own-f functional on a neighborhood
of all metrics/clocks, a controlled source-preserving local reduction,
all-order error bounds, a physical cutoff, original CD matching, scalar
health, positivity/V/G applicability, nor original P8 closure. Exceptional
zero-link or constant-root branches are treated as separate scope controls,
not silently excluded by division.
