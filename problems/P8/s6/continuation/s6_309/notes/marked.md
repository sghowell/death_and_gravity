# Marked Bose counting, total energy and independent moments

Keep exactly one marked finite residual and N additional leading-soft
emissions. For total n=N+1 identical graviton labels, the symmetric
integrand sums over the n choices of marked leg. Its current, virtual
factor and finite residual are attached to that leg's radiative state.

On the common symmetric energy simplex, relabeling maps every marked
term to the same integral. Consequently n/n!=1/N!, which is precisely
the Bose factor of the additional unmarked emissions. This remains true
when the soft kernel depends on the marked state. It does not justify
changing that state after each unmarked emission or exponentiating the
whole nonleading residual.

If the marked energy is omega and the total accepted energy is x, the
additional emissions have the total budget y=x-omega. The angularly
integrated leading measure is a_e v^(2e-1)dv. Repeated Euler beta
integration gives the N-fold simplex

[a_e Gamma(2e)y^(2e)]^N/[N! Gamma(1+2eN)].

Combining this identity with the signed-series proof gives exactly

D[R](x)=integral_(0<omega<x) P_sigma(x-omega)dR(sigma).

The signed residual has already subtracted its one-real soft singularity.
It must be integrated as that difference. Its two individual real-rate
integrals are divergent and are not separately assigned finite values.

## Independent constant-state moments

For a diagnostic kernel dR=omega^(p-1)domega, p>0, and fixed a_e,
the literal marked-energy integral at finite N is

Gamma(p)*x^p*[a_e Gamma(2e)x^(2e)]^N/
[N! Gamma(p+1+2eN)].

Summing with the same analytic virtual normalization and taking the
ordered limit gives

D_p=Gamma(p)*exp(Delta-gamma_E*a)*x^(p+a)/Gamma(p+1+a).

After division by P(x), this is x^p*Beta(p,1+a). In particular the
ratio to the unweighted marked integral x^p/p is

Gamma(p+1)*Gamma(1+a)/Gamma(p+1+a).

For positive integer p it is the product of j/(j+a), j=1,...,p.
At p=1 it is 1/(1+a), and at p=2 it is 2/[(1+a)(2+a)].
These are nontrivial suppressions for a>0. Replacing the remaining
energy by x, or applying separate individual-energy cuts, misses them.

Exact state-dependent finite-label sums, literal low-dimensional simplex
quadrature and independent full-series limits calibrate these formulas.
Diagnostic power kernels are not substituted for the original physical
47-tree residual in the theorem.
