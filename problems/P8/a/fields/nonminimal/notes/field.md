# Actual nonminimal stress, positive squares and the conformal bridge

## Flat stress and the local inequality

With eta=diag(1,-1,-1,-1), the real massless classical stress is

    T_ab=phi_a phi_b-eta_ab (partial phi)^2/2
          +xi(eta_ab Box-partial_a partial_b)phi^2.

Its off-shell trace is (6xi-1)(partial phi)^2+6xi phi Box phi.
On shell, its energy along the inertial time line is

    rho=phi_t^2/2+(1/2-2xi)|grad phi|^2-2xi phi phi_tt.

For a real compact f an exact integration-by-parts identity is

    f^2 rho
     =(f phi_t)^2/2+(1/2-2xi)f^2|grad phi|^2
       +2xi[partial_t(f phi)]^2-2xi f'^2 phi^2
       -2xi partial_t(f^2 phi phi_t).                      (1)

The boundary integrates to zero. The first three coefficients are
nonnegative for 0<=xi<=1/4. Point splitting applies (1) to the smooth
symmetric difference of any Hadamard target covariance and the vacuum.
For each real differential square, the full Fourier diagonal identity,
symmetry of the difference and positivity of the target bound the
normal-ordered diagonal below by minus the reference half-frequency
integral. The half integral has factor 1/pi.

The vacuum radial covariance is hbar/(4pi^2) integral_0^infinity
k exp[-ik(t-t')] dk. The first two groups in (1) contribute
(1-2xi)k^3. Integrating the total derivative of f phi against
exp(-i alpha t) gives i alpha times its transform, so the last positive
square contributes 2xi k alpha^2. With theta=alpha+k, the inner integral is

    integral_0^theta [(1-2xi)k^3+2xi k(theta-k)^2] dk
       =theta^4(3-4xi)/12.

For real f and transform fhat(theta)=integral exp(-i theta t)f(t)dt,
the positive-half Parseval identity is
integral_0^infinity theta^4|fhat|^2=pi||f''||^2. Thus

    integral <:rho:vac> f^2
      >= -hbar(3-4xi)/(48pi^2)||f''||^2
         -2xi integral <:phi^2:vac> f'^2.                 (2)

At xi=1/6 the coefficient is 7hbar/(144pi^2).
This agrees with the massless specialization of Theorem 4.3 of
[Fewster and Osterbrink](https://arxiv.org/pdf/0708.2450);
the local positive-squares calculation here also applies when the
target is defined only on a smooth Minkowski strip containing the
compact sample. No target-state extension to all Minkowski space
is assumed. The reference vacuum is restricted to that strip.

## Why a state qualification cannot be deleted

Fix tau>0 and an amplitude A. In a neighborhood of the sampled line use

    phi=A[3+(t/tau)^2+|X|^2/(3tau^2)].

It solves Box phi=0. At xi=1/6 the on-shell trace vanishes, and on X=0

    E=rho=A^2/tau^2[(4/3)(t/tau)^2-2]
           <=-2A^2/(3tau^2), |t|<=tau.                  (3)

The polynomial alone has infinite energy and is NOT the state witness.
Instead set its t=0 Cauchy data equal to these polynomial data on a
ball of radius 2tau, smoothly cut them to zero outside radius 3tau,
and solve the actual global wave equation. Finite propagation speed
makes the solution agree with the polynomial on a neighborhood of
the entire line |t|<=tau. The resulting compact smooth Cauchy data
have finite classical energy and finite vacuum one-particle norm:
their Fourier transforms decay rapidly, and the possible k^-1
weight is integrable with the three-dimensional radial measure.

A Weyl displacement of the vacuum by this real finite-norm solution
is a positive Hadamard coherent state. Its two-point function is
W_vac+phi(X)phi(Y), with smooth difference, and its normal-ordered
stress equals the classical stress of that solution. For every fixed
real nonzero f supported in (-tau,tau), (3) sends its averaged E
to minus infinity as A grows. Each member has finite energy; no
common finite energy bound is asserted. This proves absence of an
all-Hadamard state-independent lower bound for that average.
It neither solves SEE nor excludes state-dependent or energy-capped bounds.

## Physical conformal transport of BOTH costs

On g=a(eta)^2 eta set Phi=a^-1 phi, dt=a d eta.
The covariant improved stress is

    T_ab=Phi_a Phi_b-g_ab(partial Phi)^2/2
          +xi(g_ab Box-nabla_a nabla_b-G_ab)Phi^2.

At xi=1/6 every one of its 16 covariant components transforms as
T_g,ab=a^-2 T_flat,ab, off shell. The field equation transforms as
(Box_g+R/6)Phi=a^-3 Box_flat phi. Both are checked using arbitrary
field jets and arbitrary positive a with arbitrary first two jets.
For differences of Hadamard stresses the local anomaly cancels;
their physical E difference equals their rho difference and has a^-4
weight. The relative Wick square has a^-2 weight.

For a proper sampler f(t), put F(eta)=a^-3/2 f(t(eta)).
Direct differentiation and the measure give

    integral |F''(eta)|^2 d eta=integral |L_H f(t)|^2 dt,
    integral a^2 |F'(eta)|^2 d eta=integral |G_H f(t)|^2 dt.

The a^2 in the second identity is essential: it comes from converting
the flat relative Wick square to the physical relative Wick square.
Inserting (2) and restoring the actual reference E yields (1) of the
formulation. A one-sided upper bound w<=Phi_*^2 suffices because
the squared weight is nonnegative and its coefficient is negative.

The conformal vacuum exists by restriction and conformal transport
on any smooth positive strip. This preserves positivity, commutators
and the Hadamard condition. Arbitrary target Hadamard states remain
admitted; their smooth relative Wick square and stress on the compact
sample make the bounds continuous in H2. C-infinity compact density
therefore extends the inequality to H2_0 with both zero boundary traces.
No global upper bound on a or its derivatives is needed outside
the compact segment used for a particular inequality.
