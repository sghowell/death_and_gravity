# Complete positive soft sum and controlled regulator removal

Let x=E/nu in(0,1],e>0 and a_e>=0. After angle integration, one
leading real graviton has radial measure
a_e*nu^(-2e)*omega^(2e-1)domega.
The N-fold TOTAL-energy simplex, including its Bose factor, is

T_N = [a_e*Gamma(2e)*x^(2e)]^N /
      [N!*Gamma(1+2e*N)],     T_0=1.

This follows recursively from the Euler beta integral: integrate the
last energy between0 and E and apply the same formula at E-omega.
No upper integration bound can be dropped if it changes the total cut.

Set lambda=a_e*Gamma(2e)*x^(2e). For an AUXILIARY Poisson variable
N with mean lambda, the exact analytic-virtual-normalized sum is

exp[-a/(2e)] sum_N T_N
= exp[lambda-a/(2e)] E_Poisson[1/Gamma(1+2e*N)].

This is an identity of positive convergent sums. The auxiliary Poisson
law is a computational reweighting, not the actual conditioned
multiplicity law of the detector.

At fixed x>0 assume the established expansion
a_e=a+2e*Delta+O(e^2). The full Gamma expansion gives

lambda-a/(2e) -> Delta-EulerGamma*a+a*lnx,
E[2e*N] -> a,    Var[2e*N]=4e^2*lambda ->0.

Hence2e*N->a in probability (indeed L2). The function
g(y)=1/Gamma(1+y) is continuous fory>=0 and globally bounded:
Gamma(1+y)=integral_0^infinity exp(-t)t^y dt
>=integral_1^infinity exp(-t)dt=exp(-1).
The bounded-continuous expectation theorem, or a direct Chebyshev split
into a neighborhood of a and its complement, therefore gives
E[g(2e*N)]->g(a). This proves the complete limit

I_lead(E)/|A_hard,nu|^2
= exp(Delta)*exp(-EulerGamma*a)/Gamma(1+a)*x^a.

In particular the limit was NOT taken separately at each fixed N:
the typical auxiliary N grows as a/(2e). At a=0 the same proof holds;
at physical zero-radiation forward/backward endpoints S296 has K_e=0
for all e, so a=Delta=0 and the soft factor is one. This does not remove
the elastic Coulomb pole or its amplitude phase.

There is also a genuine positive-probability interpretation. On
0<omega<=nu use the finite-e Poisson intensity above and the normalized
virtual factor exp[-a_e/(2e)]. For E<=nu its total-energy cumulative
probability is exactly the same simplex series with that virtual factor.
The analytic virtual normalization differs by
exp[(a_e-a)/(2e)]->expDelta. In the four-dimensional limit, the normalized
energy sum has Laplace transform
exp[-a*(EulerGamma+ln(z*nu)+E1(z*nu))].
A separate numerical inversion agrees with F(a)*x^a below nu and
explicitly disagrees with an extrapolated power law above nu.

The probability construction controls the leading model. It neither
proves all-N full-amplitude soft factorization with a uniform error nor
constructs an analytic unitary scattering amplitude for the complete
physical inclusive rate.
