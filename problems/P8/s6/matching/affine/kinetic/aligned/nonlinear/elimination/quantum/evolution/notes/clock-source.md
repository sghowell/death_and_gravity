# Homogeneous clock source and physical first-derivative bounds

Use the unchanged S6.53 subtraction prescription and the
[C1 finite-integral proof](derivative.md). All variations below
have compact support inside the time interval, with the prepared
initial Gaussian data held fixed. There is no variation of the
state preparation prescription with the light background and no
claim about endpoint or infinite-time boundary terms.

## Separate the ordinary and clock-mass parts

On the actual clock a_m=b_m=1 at every time. The vector mode
equations are therefore the ordinary Proca equations, although
the physical lapse variation still contains alpha=a_N and beta=b_N.
Its pressure has the ordinary weights. Write

    rho=rho_ord+rho_mass, p=p_ord.

The ordinary energy and pressure obey the mode conservation identity.
The code checks its kinetic and potential coefficients separately
even for formal spatial dimension D. For example, with lambda=-Hz,
the cancellation conditions are -2d+D H*p_A=0 and
2lambda+2d+D H*p_B=0, with the actual transverse and longitudinal
rates and pressure weights of S6.52.

The full ordinary adiabatic coefficients through orders 0,2,4
obey this identity separately in the physical dimension. If F_n
denotes the order-2n energy bracket and P_n the corresponding
pressure bracket, the exact checks are

    D_u F_n+(1-2n)lambda F_n+3H P_n=0.

Here D_u includes the fixed-comoving momentum chain rule, not a
partial time derivative at fixed z. The S6.53 ordinary matched
local coefficients also satisfy their Ward identity. The uniform
integrable derivative bound now permits these identities to pass
through the physical subtracted momentum integral. Consequently
rho_ord'+3H(rho_ord+p_ord)=0.

## Direct clock variation, not imposed separate conservation

In normalized units the physical clock is phi=tau*u; the source
signature has x=-u'^2/N^2, with clock u'=N=1 and x=-1.
Since each mass function equals one identically on this clock,
its explicit clock-coordinate derivative at fixed x vanishes there.
Its lapse jet implies a_x=alpha/2 and b_x=beta/2. For a
homogeneous clock variation at fixed physical metric,

    delta x=-2 delta u',
    delta a_m=-alpha delta u',
    delta b_m=-beta delta u'.

The unit timelike normal is unchanged by a homogeneous monotone
clock reparametrization. Therefore the curvature coefficient tensors
in S6.53 contribute no additional first clock variation when
multiplied by their vanishing mass deviations. Terms with higher
powers of those deviations also have zero first variation.
The original source-aligned translation does not add a Gaussian
first variation: its background and first source variation vanish
as established in S6.42/50.

Let La,Lb denote the on-clock vector mass-insertion readouts per
physical volume. The mass lapse energy is
rho_mass=-alpha*La-beta*Lb. The clock variation of the same action
density is a^3*rho_mass*delta u'. Integrating by parts gives its
Euler source

    J_clock=-(rho_mass'+3H rho_mass)
           =-[rho'+3H(rho+p)].

The sign convention is delta Gamma=integral a^3 J_clock delta u
in normalized units. The code independently checks the lapse and
clock variations using time-dependent alpha,beta,La,Lb and arbitrary
higher mass powers. It does not set the mass functions to their
clock values before variation. In particular, imposing separate
ordinary-Proca conservation on rho_mass would incorrectly erase
a generally nonzero clock source.

The argument applies also to the stated counterterm density before
taking its regulator limit. S6.53 proves convergence of the finite
subtracted integrals uniformly on the compact time interval near
D=3. For compact test functions, integration by parts transfers
the source derivative to the test function, so this convergence
also justifies the limit of the clock source as a distribution.
The new physical C1 result identifies that limit with the continuous
function above. A separate uniform complex-D derivative bound is
not silently assumed.

## Local source and scale estimates

At mu=m0, the local coefficients of -J_clock relative to
1/(64pi^2), with respective m^4,m^2,1 factors, are

    -44u/[9(1+u^2)^4],
    128u(79u^2+36)/[243(1+u^2)^6],
    -64u(449u^4+898u^2+81)/[81(1+u^2)^8].

Direct differentiation of the actual S6.53 coefficients gives
these expressions. Subtracting the ordinary local energy first
gives exactly the same mass-energy balance. Thus the ordinary
vacuum energy cancels from this clock-source calculation; its
absolute value is not counted as a spurious time-dependent source.

Continuous rational envelopes on the entire I are:

| order | energy derivative | pressure derivative | clock source |
|---|---:|---:|---:|
| zero | 22/9 | 0 | 22/9 |
| two | 468535/1944 | 340/3 | 3568/243 |
| four | 447955/432 | 536 | 1186/9 |

These are exact polynomial-box bounds, with reconstruction checks,
not sampled extrema. For any of these triples U_n, the local
physical bound is sum U_n R^(4-2n)/(576 L^2), where
L=M*tau and R=m0*tau. The reference density rate and scalar-equation
scale are both M^2/tau^3. The physical clock remains phi=tau*u,
so the source has density/time units; no Planck-rescaled clock
is substituted.

Add the state derivative bound [40K/R+D1/(72R^2)]/L^2,
K=676,000,000 and D1=100,000,000, for energy and pressure.
For the clock source, add the same state energy derivative plus
6 times the sum of the S6.51 absolute state energy and pressure
bounds, since 3|H|<=6. Use the local source envelope directly,
rather than the sum of separate local absolute energy and pressure
bounds. This retains the ordinary cancellation.

At L=10^12,R=1000 the total bounds are respectively

    energy derivative:
      1913098566532399/447897600000000000000000000000,
    pressure derivative:
      5883140501/216000000000000000000000000,
    clock source:
      1494611845133383/349920000000000000000000000000.

Each is below 10^-14 in its stated reference units. The report
retains the exact local, state-evolution and reference-tail error
decomposition.

This bounds one vector Gaussian contribution at the fixed classical
clock. It does not solve the corrected light equations or control
their quantum principal symbol. Finite C1 stress/source values do
not imply all-order Hadamard admissibility, higher time derivatives,
quadratic functional variations, other field/higher loops, omitted
UV operators, an interacting cutoff, vacuum/Regge matching or V/G/B.
Original P8 is not closed.
