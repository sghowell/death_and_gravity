# Explicit cutoff-summed all-order Cauchy state

This is a new named state, not a reinterpretation of the frozen
S6.50 fourth-order Gaussian. The vector action, physical metric,
clock and S6.53 local matching prescription are unchanged.
Fix m>=1000, I=[-1/2,1/2], Amax=25/16 and
nu^2=m^2+k^2/Amax^2. Nu is time independent and equals the
physical omega at the initial slice u0=-1/2.

## All-order recurrence

Write S=1+sum_{n>=1} P_n t^n, t=omega^-2, W=omega*S,
lambda=-Hz, and D0=partial_u-2Hz(1-z)partial_z.
The logarithmic frequency rate and exact reference residual are

    L=W'/W=lambda+D S/S,
    rho_res=(1-S^2)/t-U-(D L)/2+L^2/4,
    D=D0-2lambda*t*partial_t.

For a proposed P_n, the coefficient of t^(n-1) in this residual
contains -2P_n and otherwise only earlier P_j. Indeed the new
coefficient first enters D S/S at order t^n, whereas (1-S^2)/t
contains -2P_n*t^(n-1). Thus residual cancellation uniquely
determines P_n as half the lower-order coefficient. Formal inverse
series and coefficient convolutions implement this recurrence
without treating a truncated series as convergent.

For the actual H=4u/(1+u^2), each P_n and its D0 derivative belong
to the ring of rational polynomials in u,z with denominators a
constant times a power of 1+u^2. This follows inductively because
D0 preserves that ring and the inverse formal series has leading
coefficient one. Hence every coefficient and slope has a finite,
computable rational polynomial-box envelope on I and 0<=z<=1.
The recurrence exactly reproduces the frozen P2,P4 (here P_1,P_2)
and independently cancels the leading frozen residual with P_3.

For n>=3 define

    V_n=bound(abs(P_n)),
    G_n=bound(abs(D0 P_n+(1-2n)lambda P_n)),
    C_n=max(V_n,G_n).

The physical frequency term is P_n*omega^(1-2n), so these
bounds control its value and first time derivative at once.
In particular, the sixth-order envelopes are

| polarization | V_3 | G_3 |
|---|---:|---:|
| transverse | 25941283/64 | 1424157971/64 |
| longitudinal | 12497233/32 | 705628143/32 |

## An explicit locally finite cutoff sum

Let chi(s) be zero for s<=1, one for s>=2, and between them

    chi(s)=exp[-1/(s-1)]/
            (exp[-1/(s-1)]+exp[-1/(2-s)]).

It is smooth, real and between zero and one. Smoothness follows
from flatness of b(x)=exp(-1/x) extended by zero for x<=0:
b^(j)(x)=exp(-1/x) Q_j(1/x), with polynomial recurrence
Q_0=1, Q_{j+1}(y)=y^2[Q_j(y)-Q_j'(y)].
Every such derivative tends to zero at the endpoint.

For each polarization choose

    Lambda_3=2m,
    Lambda_n=max(2Lambda_(n-1), 2^n[1+C_n]), n>=4.

Define the local frequency used solely for Cauchy preparation by

    W_B(u,k)=W_4(u,k)
          +sum_{n>=3} chi(nu/Lambda_n) P_n(u,z) omega^(1-2n).

The cutoffs depend on the fixed comoving momentum and mass, not
time, so they have zero u derivative. For each finite momentum
the sum contains only finitely many nonzero terms: thresholds at
least double. The first threshold above nu terminates the exact
preparation algorithm. Each coefficient calculation is finite.
No fixed highest adiabatic order or numerical truncation replaces
this definition.

On the support of the nth cutoff, omega>=nu>Lambda_n.
For n>=4, Lambda_n>=1 and Lambda_n>=2^n C_n imply

    C_n/omega^(2n-6) <= 2^-n.

Thus the complete frequency and slope corrections obey

    |W_B-W_4| <= F/omega^5, F=V_3+1,
    |W_B'-W_4'| <= G/omega^5, G=G_3+1.

The geometric sum of higher contributions is below one; the
actual sixth-order term is bounded separately. Since F/m^6<1/4
and W_4>=omega/2, the new frequency is everywhere positive,
W_B>=omega/4. At u0 define the exact normalized Cauchy data

    v_B=(2W_B)^(-1/2),
    v_B'=(-W_B'/(2W_B)-i W_B)*v_B,

then evolve with the unchanged exact vector mode equation.
There is no assertion that the cutoff-summed W_B solves that ODE.

For every fixed N>=3, all cutoffs up to N equal one for sufficiently
large k. The first omitted term is bounded by
C_(N+1)*omega^(-2N-1); for n>=N+2 the threshold rule bounds
each remaining term by 2^-n*omega^(-2N-1). The same argument
holds for the slope. Consequently W_B and W_B' have the complete
all-order WKB asymptotic expansion. Local finiteness gives smooth
Cauchy functions; it is not an assumption that the original formal
series converges.

## Initial-state comparison

Let W,J=W_4,W_4'/(2W_4) and V,K=W_B,W_B'/(2W_B) at u0.
The exact Bogoliubov coefficients relative to the old preparation
are

    A0=[W+V-i(K-J)]/(2sqrt(WV)),
    B0=[W-V+i(K-J)]/(2sqrt(WV)).

The value, derivative and CCR |A0|^2-|B0|^2=1 are independently
checked. Since |K-J|<=2(G+4F)/omega^6,

    |B0| <= [2F+4(G+4F)/m]/nu^6
           <= C_H/nu^6, C_H=10,000,000.

Both actual rational constants are below one million; C_H is a
declared conservative common bound. In particular |B0|<1.
The numerical error transfer uses this explicit comparison with
the frozen state, not unknown constants from a Hadamard existence
theorem. The [Proca comparison proof](hadamard.md) supplies the
all-order state property; [transfer bounds](transfer.md) quantify
the retained physical observables.
