# Full-D proper, external and endpoint calculation

## Six pairs and four external factors

S6.292 independently contracts the complete pair stress tensors in
arbitrary D and in literal D4/5/6 components. For all outgoing
on-shell momenta k_i,k_j, z=k_i.k_j and d_i=l^2+2k_i.l,
d_j=l^2-2k_j.l, it gives
N_D=4z^2-4mu^2/(D-2)+2z*l^2-2z(d_i+d_j).
The same tensor identity applies to any pair adjacent to the
nonderivative quartic. It is not a soft approximation.

The complete graph is
-C[N0_D*C0_pair+2z*B0_pair-2z*B0_on_i-2z*B0_on_j]
 /(16pi^2 kappa).
Each metric-quartic contact contributes
+C[4mu*B0_on/(D-2)+A0]/(16pi^2 kappa).
The full external factors use Z=1/(1+Sigma'), with
A0=mu(1+2e)B0_on/(1+e),
Sigma'=mu(3+2e)B0_on/[(1+e)16pi^2 kappa].

Every channel a occurs twice with2z=a-2mu. Momentum conservation
gives sum_(j!=i)z_ij=-mu. The complete sum is
-C/(8pi^2 kappa)*
 [sum_a(V_D*C0_a+(a-2mu)*B0_a)
  +mu(1+2e)*B0_on/(1+e)].
The whole radial triangle and bubble share A_a=mu-a*x(1-x).
C0_a=Gamma(-e)(4pi nu^2)^(-e)J_a(e)/2,
B0_a=Gamma(-e)(4pi nu^2)^(-e)M_a(e),
B0_on=Gamma(-e)(4pi nu^2)^(-e)mu^e/(1+2e).
This gives precisely-2T_e in the formulation.

The UV contributions of pairs, contacts and external normalization
are proportional to+4mu,-12mu,+8mu and cancel. Their separate IR
parts instead yield the full massive four-leg Bsoft. The raw
combined on-shell self-energy pole must not be read as the IR pole.

## Complete matter-metric endpoints

In each endpoint bubble retain its entire stress numerator before
the D limit. With Delta=mu-a*x(1-x), isotropic loop averaging gives
an eta coefficient
(2/D-1)(A0+Delta*B0)+(mu+a*x(1-x))*B0.
Using the exact relation A0=2Delta*B0/(D-2), this equals
2a*x(1-x)*B0. The qq coefficient is-2x(1-x)*B0.
Thus the full bubble is-2x(1-x)*Q*B0, Q=qq-eta*a.
The metric-quartic matter tadpole and full covariant on-shell mass
countervertex cancel as a whole, leaving this transverse tensor.

The literal factors iC, -iT, i^2, iB0/(16pi^2) and the identical
internal-matter half give
deltaGamma=-C*Q*integral x(1-x)B0/(16pi^2).
For the opposite on-shell tree endpoint,
T=2PP-Q/2, P.q=0, P^2=mu-a/4.
Its trace is2mu+(D-2)a/2. The full harmonic contraction is
Q.P.T=a*trace(T)/(D-2),
so two endpoint assignments produce
2Q.P.T/a=a+4mu/(D-2).
Combining every crossed channel yields the stated+E_e term.

No pure pole projection or four-dimensional-only trace was used.
The complete scalar cut in notes/cuts.md independently checks the
whole nonlocal proper-plus-endpoint amplitude and its normalization.

## Local and literature boundaries

E0=sum_a(a+2mu)/6=5mu/3 because s+t+u=4mu.
The raw local UV pole is therefore-5Cmu/(48pi^2 kappa e) in the
retained on-shell/curvature organization. Finite local matching is
not thereby chosen. E1 includes the full evanescent local term-mu
in addition to sum_a(a+2mu)integral y ln A_a; this constant has zero
b20 but is not removed from the full amplitude.

The original PDF of Gonzalez-Martin and Martin,
https://arxiv.org/pdf/1711.08009, printed page8 equation3.6 was
visually inspected. Its displayed1/12*[1/2(s+t+u)+m^2] evaluates
to m^2/4 on shell, whereas the next displayed line uses m^2/2.
That internal discrepancy prevents using the coefficient as an
independent normalization benchmark here. We assert no correction
to the paper and no agreement with its numerical UV coefficient.
The present tensor, on-shell and complete-cut derivations stand
independently; a physical beta-function conclusion is not claimed.
