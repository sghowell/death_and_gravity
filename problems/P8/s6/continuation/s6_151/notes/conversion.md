# Complete family conversion at a common regulator

For the core, let I(P) be the outer bubble and I0=I(0).
The old four-term forest is

    Old(P)=J_bare(P)-I0 I(P)-J_bare(0)+I0^2.

Let K be the momentum-independent overall pole part after
proper MS subtraction. Then

    New(P)=J_bare(P)-I_MS I(P)-K,
    New(0)=J_bare(0)-I_MS I0-K.

Literal subtraction before any limit gives

    New(P)-Old(P)=New(0)+(I0-I_MS)[I(P)-I0].

The last bracket is already finite; hence in the limit this
is F+(ell/Q) I_R(P). Evanescent pieces matter in F, but cannot
multiply a pole in the finite last bracket. This distinction
prevents either premature truncation or spurious finite terms.

The sixteen full vertex numerators sum to (C(P)+T_A)C(Q)^2,
where T_A contains BOTH internal-heavy exchange terms.
Split this, as in S6.124, into C(P)L^2 and

    H(q)=L^2 T_A+[C(P)+T_A][C(Q)^2-L^2].

The latter has a decaying full heavy propagator and a convergent
outer integral. Changing its proper inner bubble from I-I0 to
I-I_MS therefore adds precisely (ell/Q) integral H/(D1 D2).
The complete family conversion, with total assignment weight
three, is consequently the old family plus

    L^2 F sum_channels C(P)
    +(ell/Q) sum_channels {
       L^2 C(P) I_R(P)+integral H(q)/(D1 D2) }.

Both heavy triangles, the full heavy denominators and all
external routings are retained. An anchor-only conversion
would miss the second line.

This is an assigned interaction-forest conversion at fixed
canonical reference coordinates. It must not be counted again
as an independent parameter-map contribution. The remaining
global scale/coupling/field re-expansion is still a separate
assembly problem.
