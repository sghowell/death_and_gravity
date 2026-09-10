# Actual overlap and local propagator disc

Use precisely the actual S6.110 parameters and S6.111 coefficient

    lambda=10^-600, gamma=1024 times 10^-800,
    c=4 lambda^2/gamma.

The numerator N=c/2-2lambda+9gamma/4 is strictly positive.
Since pi>3, the one-loop on-shell mixing satisfies

    0 < w=N/(16pi^2) < N/144 < 10^-404.

Thus the mapped one-loop residue is 1-2w>1-2 times 10^-404>0.
The decimal value about 1.24 times 10^-405 is a diagnostic only;
the certificate uses rational inequalities.

On |s-1|<=1, B=-s has |B|<=2. The finite composite polynomial's
coefficient triangle inequality gives the simultaneous bound

    |W_MS(-s)| < U=(c/2+5lambda+5gamma)/144 < 10^-404.

The parent has D_Phi(s)=s-1+hbar Pi_R(s), with
|Pi_R(s)|<=epsilon_parent |s-1|^2 and
epsilon_parent<10^-405. The inverse-map two-point function,
through one loop, is

    G_Psi = 1/(s-1) - hbar Pi_R(s)/(s-1)^2
            - 2 hbar W_MS(-s)/(s-1).

Equivalently its inverse at the same order is

    D_Psi = s-1 + hbar[Pi_R(s)+2(s-1)W_MS(-s)].

Multiplication and the pole value/derivative are checked exactly.
For the displayed one-loop inverse at hbar=1, factoring s-1
leaves a bracket whose distance from one on the unit disc is
at most epsilon_parent+2U<3 times 10^-404<1. Hence it has no
additional zeros there. The derivative gives residue 1-2w
through one loop, consistently with composite pole factorization.

All statements use the same local disc and fixed subtraction
prescription. They neither control the exact higher-loop
propagator nor classify zeros outside this disc. Treating this
polynomial correction as a new exact free higher-derivative
kinetic theory would not be the perturbative model used here.
