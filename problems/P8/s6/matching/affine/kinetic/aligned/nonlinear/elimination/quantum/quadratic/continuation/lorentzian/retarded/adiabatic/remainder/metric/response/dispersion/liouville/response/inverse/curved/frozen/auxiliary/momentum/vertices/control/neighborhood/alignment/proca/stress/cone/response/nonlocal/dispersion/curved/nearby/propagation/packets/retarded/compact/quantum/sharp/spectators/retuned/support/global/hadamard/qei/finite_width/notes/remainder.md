# Whole-interval graph and exact-mode error

For every u in [-h,h], h=1/100, the ball-jet variable is u+z,
where the constant coefficient encloses the ENTIRE real interval
and the first z coefficient is one. Polynomial convolution, pivoted
reciprocal and the positive square-root recurrence enclose the
actual Taylor coefficients at every such u. A failed pivot or root
domain is rejected. Differentiation shifts factorial-normalized
coefficients. This is interval automatic differentiation, not a
finite collection of time samples.

Order-six Riccati coefficients are computed with nine time-jet
orders, retaining three derivatives of every R_n. They use the
literal h0,h1,h2 and the same unique frequency-sum inverse as the
exact field calculation. Symbolic induction cancels orders zero
through six; interval bounds are used for the NONZERO remaining
coefficients, never to certify a zero from an interval containing it.

For R6_trunc, write its scaled Riccati defect as
sum_(n=7)^14 d_n k^(-n). The explicit coefficient enclosures give,
for k>=K=32,

    ||defect|| <= 91215897 k^(-7).

Put x=V z. The approximate configuration equation is exactly

    z'=(-ik omega+B(u,k))z,  B=sum_(n=0)^7 B_n k^(-n).

All coefficient derivatives through order three are bounded.
Their finite weighted sums at k=32 give

    ||B^(j)|| <= (2,24,300,7082)_j.

The positive diagonal frequencies are >0.99, with derivative
upper bounds (2,1,5,103). In particular
||z(t)||<=exp(2h)||z(0)|| on the initial strip; the principal
oscillation makes no contribution to the real norm derivative.

Factor sqrt(hbar/(2 kappa))*sqrt(k) out of physical packet modes.
For each normalized initial column, approximate with initial
configuration B_B(0,k)^(-1/2) but initial momentum
R6_trunc(0,k) B_B(0,k)^(-1/2). Thus the nonzero initial error is
retained and bounded by 1200002 k^(-8). The approximate mode
configuration has norm <=4 exp(2h), while its full phase residual
is purely in the momentum block, with norm
<=4 exp(2h)*91215897 k^(-6).

Duhamel with the COMPLETE original-chart transfer, using the safe
duration 2h and sqrt(182)<14, proves

    ||Y_exact-Y_app|| <= C_E k^(-6),
    C_E=14 exp(196/25) [600001/512+91215897 exp(1/50)/25].

The configuration-row projection of the phase residual vanishes.
Therefore the physical time derivative error is bounded by
16 C_E k^(-11/2), and the spatial-gradient sum by
C_E k^(-11/2), after removing sqrt(hbar/(2 kappa)).
This step retains all time-dependent row factors; it does not
differentiate an O(k^(-6)) statement without its evolution equation.
