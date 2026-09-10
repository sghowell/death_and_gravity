# Proper and overall subtractions at the same regulator

Let I(Q) denote the regulated inner light bubble and
I0bar=I0/(16pi^2) the entire same mass-one,
zero-external-momentum reference of S6.113.
The proper subtraction gives

    I_R(Q)=I(Q)-I0bar
          =-(16pi^2)^-1 integral_0^1
                         ln[1+a Q_E^2] dx, a=x(1-x).

The local proper counterterm is defined before
attaching the heavy propagators adjacent to that
core. Its external-field type is taken from the
actual full two-field graph.

For each of the two overall cases let J_bare(P)
be the unrenormalized outer integral, J_bare(0)
its zero-external-momentum value, and I(P) the
outer light bubble. Its four restricted forests give

    J_bare(P) - I0bar I(P) - J_bare(0) + I0bar^2.

The final term is the nested subtraction applied
inside out. Define at this same regulator

    J_R(P)=J_bare(P)-I0bar I(P),
    J0=J_bare(0)-I0bar^2.

The result is J_R(P)-J0. This explicitly fixes the
new-order recursive zero-momentum reference.
It is a prescription, not something determined
by the one-loop conditions alone. J0 is not
separately assigned a finite number.

The other fourteen forests perform only the proper
subtraction; their remaining integrals are finite.

At an outer channel write

    C(P)=-L+g/(M-z),  V_A=C(P)+T_A,
    C(Q)=-L+g/(M+Q_E^2).

Here T_A is the sum of the two full internal-heavy
exchange terms, not its momentum expansion.
The 16 selected numerators sum exactly to
V_A C(Q)^2. Split this as

    C(P)L^2 + L^2 T_A + V_A[C(Q)^2-L^2].

The renormalized family is consequently the sum of

    C(P)L^2 [J_R(P)-J0]

and the outer integral of

    [L^2 T_A+V_A(C(Q)^2-L^2)] I_R(Q)/(D1 D2).

The identities are applied at a common regulator
before its removal, so products of divergent and
evanescent pieces are never discarded independently.

The overall subtraction is -C(P)L^2 J0 per channel.
Across all channels it is implemented by the
literal local full-model tree variations

    delta L=-3L^3 J0,
    delta g=-L^2g J0,
    delta M=0.

These are this graph family's counterterm contributions.
The loop measures are already included in J0.
Other fixed finite-potential insertions and complete
two-loop normalization are not inferred from them.
