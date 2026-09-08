# Exact subtraction tail and finite physical mode integrals

Retain S6.50, its unchanged action, physical lapse/pressure forms and
exact Gaussian initial data. All bounds below hold on I=[-1/2,1/2]
for every comoving momentum and m=m0*tau>=1000. Work in tau=1 units
until the final scale restoration. No previously chosen finite
counterterm is overwritten.

## Physical reference and derivative-order subtraction

Use omega^2=m^2+k_com^2/a_s^2, z=(omega^2-m^2)/omega^2 and
t=omega^-2. In this proof t is not time. Write p=P2, r=P4,
b2=B2, b4=B4 and S=1+pt+rt^2 for the S6.50 reference.
Both energy and pressure per polarization can be written

    E(v)=[A|v'-dv|^2+B*omega^2*|v|^2]/(2a_s^3).

The physical coefficients, with alpha=4/(9h), beta=28/(81h), are

| observable | A | B | c1=d+lambda/2 |
|---|---|---|---|
| transverse energy | 1 | 1+beta(1-z) | H(1-z)/2 |
| longitudinal energy | 1-alpha*z | 1+beta | H(1+z)/2 |
| transverse pressure | 1/3 | (2z-1)/3 | H(1-z)/2 |
| longitudinal pressure | (1+2z)/3 | -1/3 | H(1+z)/2 |

Here lambda=omega'/omega=-Hz and h=(1+u^2)^3. Four independent
symbolic comparisons recover the original S6.50 quadratic forms,
including the actual mass jets. In particular the pressure coefficients
are already divided by three; no extra polarization or pressure factor
may be inserted in the following common normalization.

For the normalized reference f, E(f)=omega*F/(4a_s^3), exactly, with

    F=A*S+B/S+A*t*[c1+(b2*t+b4*t^2)/(2S)]^2/S.

To assign adiabatic order introduce a separate formal epsilon:
p->epsilon^2 p, r->epsilon^4 r, b2->epsilon^3 b2,
b4->epsilon^5 b4, c1->epsilon c1. The mass, omega, z, A and B
have order zero. This is derivative counting, not a large-mass
expansion at incorrectly fixed physical momentum. Expansion through
epsilon^4, independently checked by symbolic series, gives

    F0=A+B,
    F2=(A-B)p+A*c1^2,
    F4=(A-B)r+B*p^2+A*(c1*b2-c1^2*p),
    E_AD[0:4]=omega*(F0+t*F2+t^2*F4)/(4a_s^3).

The complete coefficients are subtracted, including any parts whose
momentum integrals are finite. No unmentioned finite local remainder
is added. The terminology and the ordinary-Proca cross-check below
follow the derivative-order prescription discussed in
[Maranon-Gonzalez and Navarro-Salas](https://arxiv.org/html/2412.01963v1).
Their result is not imported as a covariant matching theorem for the
present clock-dependent mass tensor.

## Exact rational tail

Define the following rational functions:

    R1=[2pr-p^3+r(r-p^2)t]/S,
    R2=[p^2-r+prt]/S,
    R3=[b4-b2(2p+(p^2+2r)t+2prt^2+r^2t^3)]/S^2.

Direct multiplication and differentiation-order expansion give

    F-F0-tF2-t^2F4=t^3*R,
    R=B*R1+A*c1^2*R2+A*c1*R3+A*(b2+b4t)^2/(4S^3).

This is an identity at every allowed frequency, not an asymptotic
remainder claim. S6.50 proves S>1/2 and supplies exact continuous
coefficient bounds P>=|p|, Q>=|r|, B2>=|b2|, B4>=|b4|.
Set tmax=10^-6. On the full domain, |A|<=1, |B|<3/2 and
|c1|<=2. These give

    R1_abs=2*[2PQ+P^3+Q(Q+P^2)tmax],
    R2_abs=2*[P^2+Q+PQtmax],
    R3_abs=4*[B4+B2(2P+(P^2+2Q)tmax
                              +2PQtmax^2+Q^2tmax^3)],
    |R|<=(3/2)R1_abs+4R2_abs+2R3_abs+2(B2+B4tmax)^2
        <D, D=400000.

Both polarization inequalities use the frozen rational box envelopes
and are checked exactly. Therefore the same envelope, including signed
pressure, holds per polarization:

    abs(E(f)-E_AD[0:4])<=D/(4a_s^3*omega^5).

There is no ultraviolet or infrared cutoff. With nonzero mass the
physical-momentum change of variables gives

    integral d^3k/(2pi)^3 omega^-5=a_s^3/(6pi^2*m^2),

using integral_0^infinity y^2/(1+y^2)^(5/2) dy=1/3, independently
evaluated. Sum two transverse and one longitudinal polarization to get

    abs(reference subtracted energy or pressure)<=D/(8pi^2*m^2)
                                                   <D/(72m^2).

This also bounds the integral of the absolute value and is uniform
in time. Thus the UV cutoff can be removed from the difference without
assigning separate finite values to any divergent term.

## Exact Gaussian observables in this stated subtraction prescription

Define E_sub(v) to be the full momentum integral of
E(v)-E_AD[0:4], with the S6.50 exact modes v. Split its integrand
as [E(v)-E(f)]+[E(f)-E_AD[0:4]]. Both brackets have uniform
absolutely integrable envelopes already proved. S6.50 has C=2000000,
so the resulting finite observables satisfy

    abs(rho_sub)<=12C/m+D/(72m^2),
    abs(p_sub)<=108C/(5m)+D/(72m^2).

Restoring units and dividing by M^2/tau^2 gives the same expressions
with m replaced by R=m0*tau and a further factor (M*tau)^-2.
At M*tau=10^12,R=1000 the reference-tail contribution is at most
1/(180*10^24), and the complete subtracted energy/pressure bounds
are 4320001/(180*10^24) and 7776001/(180*10^24), respectively.
The bound is for this explicit prescription. Arbitrary finite local
counterterms have not been bounded or set to zero in another scheme.

## Independent integrated and failing-subtraction controls

As a separate alpha=beta=0 ordinary-Proca control, integrate the two
transverse second- and fourth-order coefficients exactly. Express the
polynomial in y=1-z=m^2/omega^2. Every UV-convergent monomial has
the elementary radial Beta/Gamma integral. A term outside that
convergence domain is rejected, not analytically assigned a finite
integral. The result is

    rhoT2=m^2*H^2/(48pi^2),
    rhoT4=[H'^2-6H^2H'-H^4-2HH'']/(480pi^2).

These independently reproduce the integrated transverse terms in
Eq. (50) of the primary paper cited above. This benchmark includes
the two polarizations and our Wronskian-i normalization. It is not
the actual clock-sensitive lapse result with nonzero alpha,beta.

For that actual result, at u=0 the high-momentum limit of
2F4_transverse+F4_longitudinal is -4/27. Therefore subtracting
only orders zero and two leaves -log(Lambda)/(54pi^2) at the
bounce, plus terms with finite UV limits. S6.50's exact evolution
difference is finite and cannot cancel this divergence. The complete
fourth-order subtraction is essential for the present observable.

## Boundary of the result

This supplies finite, explicitly subtracted Gaussian energy and
pressure on the actual compact rolling interval. It does not yet
identify the full finite covariant mass/metric counterterms that
realize this prescription, nor prove all-order Hadamard admissibility
of the selected finite-order initial data. Agreement of one integrated
ordinary-Proca control or of the zero-order pole is insufficient for
that matching. No conservation identity is borrowed for a vector
component that exchanges with the clock equation.

All-time state control, other field loops, higher-order loops, the
quantum-corrected bounce and constrained cones, interacting cutoff
and original V/G/B admissibility remain unproved. Original P8 is open.
