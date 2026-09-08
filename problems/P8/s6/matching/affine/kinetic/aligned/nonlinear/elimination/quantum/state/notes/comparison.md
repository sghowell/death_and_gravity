# Exact curved modes and a finite all-momentum energy difference

Use the unchanged canonical equations of S6.42, at fixed comoving
momentum on I=[-1/2,1/2], and normalized mass m=m0*tau>=1000.
The physical scale factor a_s=(1+u^2)^2 lies in [1,Amax],
Amax=25/16, and |H|<=2. These are the actual rolling equations,
not a frozen-frequency, de Sitter or Euclidean surrogate.

Write omega^2=q+m^2, q=k_com^2/a_s^2, z=q/omega^2 and t=omega^-2.
The variable t in the symbolic program is inverse frequency squared,
not time. At fixed k_com and m the total u derivative is

    D=partial_u-2H*z*(1-z)*partial_z-2lambda*t*partial_t,
    lambda=omega'/omega=-H*z.

Independent chain-rule substitutions check this against the original
canonical frequency, including the longitudinal correction

    UT=H'/2+H^2/4,
    UL=(1/2+z)H'+(1/4-z+3z^2)H^2,
    v''+(omega^2-U)v=0.

## Positive reference and exact residual, not an asymptotic estimate

Define W=omega*S, S=1+P2*t+P4*t^2, where

    P2=-U/2-Dlambda/4+lambda^2/8,
    P4=-P2^2/2-D^2P2/4+5lambda*DP2/4
                     +(Dlambda/2-3lambda^2/2)*P2.

Introduce B2=DP2-2lambda*P2, B4=DP4-4lambda*P4,
A=P2^2+2P4 and B=-(DB4-4lambda*B4)/2+lambda*B4/2.
For f=(2W)^-1/2 exp(-i integral W du), its residual is exactly

    f''+(omega^2-U)f=rho*f,
    rho=t^2*[-A*(P2+P4*t)/S-2P2P4-P4^2*t+B/S
                +3*(B2+t*B4)^2/(4S^2)].

The program independently differentiates W'/W=lambda+DS/S and
checks this entire rational identity and both lower-order cancellations.
Discarding P4 leaves a nonzero term at order t; that approximation
does not have the residual bound used below.

Every required coefficient has the form P(u,z)/(1+u^2)^n. If
P=sum c_ij u^i z^j, then |P|/(1+u^2)^n<=sum |c_ij|/2^i
throughout the closed box |u|<=1/2, 0<=z<=1. Denominator positivity
and exact polynomial reconstruction are checked. This is a continuous
bound, including arbitrarily large momentum, not a grid test.

| coefficient | transverse envelope | longitudinal envelope |
|---|---:|---:|
| P2 | 7 | 6 |
| P4 | 13843/16 | 3219/4 |
| B2 | 119 | 104 |
| B4 | 262453/8 | 129761/4 |
| A | 6753/4 | 3249/2 |
| B | 24893281/32 | 12828497/16 |

Using t<=10^-6, these imply |S-1|<1/2 and |W'/W|<4.
Apply S^-1<=2 and S^-2<=4 in the exact residual. The two bounds
on |rho|/t^2 are below 1634089 and 1665185 respectively; hence

    omega/2<W<3omega/2, |rho|<=C/omega^4, C=2000000.

All constants and inequalities are rationally certified. They are
deliberately generous; no unproved WKB remainder theorem is used.

## Exact normalized modes and variation of constants

Prepare the exact mode v at u0=-1/2 by v=f, v'=f'. Real smooth
coefficients give a unique exact solution on I for every momentum.
The canonical Wronskian is f f*'-f* f'=i; real exact evolution keeps
it i. These modes define a positive normalized Gaussian comparison
on the physical polarization algebra. Smooth massive low-momentum
data and polynomial high-momentum bounds give distributional
two-point functions. The homogeneous vector chart handles k_com=0.

Write v=A_k f+B_k f* with A_k' f+B_k' f*=0, initially (A_k,B_k)=(1,0).
The exact first-order system is

    [A_k',B_k']^T =
      [[-i*rho*|f|^2, -i*rho*(f*)^2],
       [ i*rho*f^2,    i*rho*|f|^2]] [A_k,B_k]^T.

The program checks the readout constraint, exact mode equation and
pseudo-unitarity with metric diag(1,-1), so |A_k|^2-|B_k|^2=1.
The matrix column-sum norm is |rho|/W. Set

    nu_k^2=m^2+k_com^2/Amax^2,
    J_k=integral_I |rho|/W du <=2C/nu_k^5.

Gronwall and the integrated equation give |A_k|+|B_k|<=exp(J_k)
and |A_k-1|+|B_k|<=exp(J_k)-1. The triangle inequality applied
to v-f then gives

    abs(|v|^2-|f|^2)<=(exp(2J_k)-1)|f|^2.

The same estimate holds for v'-d*v and f'-d*f for any real d,
because the variation constraint cancels coefficient derivatives and
conjugation preserves the reference norm. Since J_k<=4*10^-9<1/4,
the positive Taylor series gives exp(2J_k)<=1/(1-2J_k), and thus

    exp(2J_k)-1<=4J_k<=8C/nu_k^5.

## Energy and pressure integrals

Use the physical forms from [energy.md](energy.md), including the
clock mass jets. They have positive A,B and |d|<=3. The exact
reference derivative norm is

    |f'-d*f|^2=[W^2+(d+W'/(2W))^2]/(2W).

The frequency-scaled energy numerator is bounded by
9/4+3/2+25/m^2<4, so rho_pol(f)<=2omega/a_s^3.
Consequently, per polarization,

    abs(rho_pol(v)-rho_pol(f))<=16C*omega/(a_s^3*nu_k^5)
                             <=16C*Amax/(a_s^3*nu_k^4).

There are precisely three polarizations. Using the full R^3 measure,

    integral d^3k/(2pi)^3 /nu_k^4 = Amax^3/(8*pi*m),

as independently checked by the elementary radial integral
integral_0^infinity y^2/(1+y^2)^2 dy=pi/4, gives uniformly on I

    abs(Delta rho)<=6C*Amax^4/(pi*a_s^3*m)<12C/m,
    abs(Delta p)<=108C/(5m).

The latter uses the absolute pressure weights, not a positivity
assumption for pressure. These are absolutely convergent integrals
of differences; neither individual unsubtracted energy is assigned
a finite value. The envelope also bounds the absolute-integral
tails and justifies passing from finite momentum cutoffs to infinity.

Restoring units multiplies normalized densities by tau^-4. Relative
to M^2/tau^2, with L=M*tau and R=m0*tau, the bounds are
12C/(R*L^2) and 108C/(5R*L^2). At L=10^12,R=1000 they are
2.4*10^-20 and 4.32*10^-20 respectively. This does not assert a
small total renormalized energy: the reference has not yet been
subtracted or matched to a covariant prescription.

## State and quantum scope

The prepared fourth-order data and comparison satisfy the displayed
CCR and energy estimates. Finite adiabatic order is not being called
an all-order Hadamard condition. The distinction matters for local
observables; see [Junker and Schrohe](https://arxiv.org/abs/math-ph/0109010).
For vector-specific adiabatic constructions and constrained propagator
issues see [Maranon-Gonzalez and Navarro-Salas](https://arxiv.org/html/2310.11860v2).
No theorem about a scalar field is silently substituted for this
constrained vector. The exact oscillator proof above is self-contained.

The reference's full adiabatic subtraction, covariant mass/metric
counterterm matching and an admissible state prescription for all
required observables remain to be supplied. This is a compact-time,
vector-only evolution difference, not an all-time state estimate,
finite total quantum stress, corrected bounce/constraints/cones,
interacting cutoff, higher-loop bound, V/G/B certificate or P8 closure.
