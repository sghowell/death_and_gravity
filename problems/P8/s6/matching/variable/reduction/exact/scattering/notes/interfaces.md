# Action, domain and independent-audit interfaces

This child uses the frozen S6.33 exact isotropic own-f branch at the same
off-shell physical metric g=eta and clock theta=T. It does not use the
different rolling background or deform the parent interactions. S6.34 is
a sibling: its shrinking-duration pulse is neither an assumption nor the
source of the present fixed-window conclusion.

## Literal action and physical units

The frozen norm-two tensor density is

    L_T = M^2 q_T^2/4 + M^2 k Q_T^2/4 - 2 beta1 b (Q-q)^2/4.

With u=T/tau and beta1=(M^2/tau^2)B1, the action is M^2/tau times
the integral of

    L_u = q_u^2/4 + k Q_u^2/4 - spring (Q-q)^2/4,
    k=b^3/N, spring=2B1 b.

`bridges.identities()` pulls the frozen tensor density back exactly and
compares it with the independently authored `action.py`. It also compares
the rational background profiles and the desingularized F and b-cubed
functions. The implicit branch and positive root are inherited from the
recursively rebuilt report; the jet derivatives are not point samples.

Let n=sqrt(k) and Y=nQ. Substitution into L_u, followed by one integration
by parts, gives

    L_normal = q_u^2/4 + Y_u^2/4 - V Y^2/4
               + spring Y q/(2n) - spring q^2/4,
    V = spring/k - n_uu/n,
    L_u[Q=Y/n] - L_normal = D_u[-n_u Y^2/(4n)].

Hence the exact canonical equation is Y_uu+VY=(spring/n)q. This source
weight is not spring/k. Under the second change

    r=sqrt(u^2+delta/8), t=asinh(sqrt(8)u/sqrt(delta)), Y=sqrt(r)psi,

the right side is multiplied by r^(3/2). `action.identities()` checks
both weights and the entire pump, even though the theorem here sets q=0.
No conclusion about a nonzero source follows solely from these identities.

The physical potential is V/tau^2. Its remainder bound is 44/tau^2,
not 44 without units. At an endpoint, Q_u=tau Q_T, so the physical-velocity
map is obtained by multiplying the dimensionless map's second column by
tau. Its determinant is tau k/rho>0. Four entrywise identities compare
the independently derived action and connection endpoint maps, retaining
the actual k_u term. These changes normalize the hidden tensor only;
the physical metric and its matter coupling are not redefined.

## The connection theorem's premise is discharged

`potential.py` proves |V-16/(delta+8u^2)|<44 on |u|<=1/100 and
0<delta<=1/100. `connection.py` proves a conditional homogeneous transfer
comparison on the same u interval but only 0<delta<=10^-6. The wrapper
replays the former proof, equates the constants, checks the domain subset,
and then replays every latter margin. Thus the combined certificate is
not conditional on an unverified coefficient assertion.

The fixed physical duration is tau/50. The coefficient-extraction frame
is unitary in (psi,psi_t/rho), and the fictitious free exterior defines
exactly the finite-endpoint matrix. It supplies no metric, source or
vacuum outside the certified interval. The analytic Duhamel estimate
includes both discarded reference tails. The resulting bound on a
matrix element is basis-specified, not an assertion that every physical
input has a particular response.

## Independent calculations and their limits

The primary engine differentiates in two independent variables using
exact Fraction intervals and total-degree-three nilpotent Taylor algebra.
`arb_audit.py` instead solves the implicit Taylor coefficients one at a
time in a univariate Arb algebra. It has its own reciprocal, cubic-root
and derivative recurrences and does not import the primary jet engine.
An exact eight-box cover of the complete (v,delta,zeta) box proves the
same |A_v|<129 and |pump|<22 bounds. This is a cover, not a sampled grid.
Independent symbolic tests check polynomial operations and exact center
derivatives at several precisions.

`connection_audit.py` adds independent Acb special-function enclosures.
Overlapping or zero-containing identity residuals are numerical
corroboration, not a proof that a residual is identically zero. The
analytic Gamma/Gauss source audit and written ODE connection prove the
identity. Strict ball inequalities can separately corroborate finite
positive margins. No approximate decimal is promoted to exact report data.

The source-pinned report recursively rebuilds S6.33 and its unchanged
ancestry. It checks exact residuals, continuous-domain proof margins,
independent enclosures and invalid-input controls. The result is a written
mathematical certificate with independent computation, not proof-assistant
formalization, a full source-matched EFT or original P8 closure.
