# General finite-difference proof of the IR-divided generating function

Writec0=ln(tau)-ln(4pi)+gamma_E andc=c0+2/V.
The exact Gamma expansion for|argument shift|<1 gives

log H_L=L*e*c0+sum_(j>=2) zeta(j)*P_j(L)*e^j/j,
P_j(L)=(L+1)(-1)^j+L^j-(-1)^j(L+1)^j.

The polynomialP_j has degree at mostj.
Its leading coefficient is2 for oddj and0 for evenj.
Also

log(eta_e/eta)
=log[1+(1+2/V)e]-log(1+e).

Multiplying this logarithm byL gives first coefficient2L/V and higher
e^j coefficients of degree1 inj's polynomial variableL.
Thus, in
f_L(e)=(eta_e/eta)^L*H_L(e),
the coefficient ofe^r is a polynomial inL of degree at mostr.
This follows directly by expanding the exponential: the sum of
degrees is bounded by the sum of the corresponding powers ofe.

At coupling ordereta^N, multiplication by the SAME IR divider
exp(-i*eta/e) gives

i^N/e^N * sum_(L=0)^N
 (-1)^(N-L)*f_L(e)/(L!*(N-L)!).

For any polynomialP of degree<N this factorial-normalized Nth
difference vanishes. For degreeN it equals the leading coefficient.
One proof applies the difference operatorN times: each application
lowers degree byone and multiplies the leading coefficient by the
previous degree. Dividing byN! gives the asserted normalization.

Every negative power ofe therefore cancels for arbitrary fixedN.
The finite term only selects the highest possibleL-degree.
In the logarithm that selection retainsLc at ordere and
2*zeta(j)*L^j/j for oddj>=3; all evenj leading terms cancel.
All higher derivatives ofV_e remain in the intermediate expression
but have degree too small to survive this specific selection.

Consequently the complete finite coefficient generating function is

exp[i*eta*c+2 sum_(odd j>=3) zeta(j)*(i*eta)^j/j]
=exp(i*eta*c)*exp(-2i*gamma_E*eta)
 *Gamma(1-i*eta)/Gamma(1+i*eta).

The equality converges for|eta|<1 at every fixedtau>0.
It proves an all-order statement about the stated leading graph class,
not about convergence of the entire interacting perturbation series.

The first coefficients are1, i*c, -c^2/2 and
-i*(c^3/6+2*zeta(3)/3). Independent exact checks retain the complete
mass-ratio and Gamma polynomials throughsix loop orders; independent
raw-Gamma numerical limits confirm these finite coefficients after
all pole cancellations. Those finite checks supplement the general
degree argument and are not its justification.

Atone loop the coefficient is
i*V/(8pi*kappa*D)[ln tau-ell+2/V], exactly the separately assembled
S303 principal phase. DroppingV_e's derivative would miss2/V.
