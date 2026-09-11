# An explicit fixed coefficient, not an adaptive subtraction

Let rho_ref(u),P_ref(u) be the actual ordinary Proca physical
stress of the S6.176 reference problem. Its definition is the
full radial exact-minus-fourth-order-adiabatic integral plus
the same finite covariant local matching terms, with all three
polarizations and exactly the S6.55 all-order Cauchy data at
t0=-1/2. The source integrand definition is retained in the
report. No finite momentum integration cutoff is introduced.

Set r=rho_ref/kappa, p=P_ref/kappa and
p_v=5m^4/(128pi^2 kappa), m=1000, kappa=1e800.
The positive-frequency Minkowski mode is exactly its
order-zero adiabatic reference; the remaining vacuum pressure
is p_v and its energy is -p_v in the fixed prescription.

Using the same integer N=1024, define the COMPLETE coefficient

    T(X)=X^N/[X^N+(1-X)^N],
    DeltaF(u,X)=-p_v+T(X)[-p(u)+p_v
                              -(r(u)+p(u))(X-1)/2].

The new candidate is CD-REG-AFFINE-ISO-QG1, with
S_new=S_old+kappa integral sqrt(-g) DeltaF.
This is a finite physical coefficient choice, NOT a change of
renormalization scheme, state-dependent normal ordering, or a
claim that an old action had this term.

The reference functions are computed/defined ONCE from the
specified history and state, then frozen. On another metric,
DeltaF still uses these same functions of its scalar argument
u; it does not recompute the stress of that metric. A function
chosen using a reference solution is a local coefficient,
not a functional of the live history or its future evolution.

The even-power denominator is strictly positive for every
real X. The reference free stress is smooth on every compact
time interval. Thus the new coefficient is globally smooth
in real u,X; real analyticity in u is not asserted.

At X=1, 1-T=(1-X)^N/[X^N+(1-X)^N]. Hence through X order1023
the coefficient agrees exactly with the affine clock profile

    A(u)+B(u)(X-1), A=-p, B=-(r+p)/2.

In particular its second X derivative vanishes on the clock,
but the whole open tube is NOT falsely identified with a
linear polynomial. Its small off-clock remainder is bounded
separately.

DeltaF is independent of the unrestricted affine connection,
W and the original M1 field. The exact60-quotient reduction,
56 algebraic complement,4 projective gauge directions, R,
the regular lower dictionary, full source S, physical metric
chart, vector mass and canonical conditional Gaussian operator
are unchanged. Only the already-enabled scalar F coefficient
is changed. No disabled function group is activated.
