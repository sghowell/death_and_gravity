# Finite current bound and labeled coefficient majorant

## Statement and domain

Use S317's complete temporal pure-soft recursion, unchanged canonical
Einstein vertices, physical future null free rays of energies w_i>0,
and unit-Frobenius spatial TT leaves. Work away from exact internal
propagator poles. For a nonempty subset S let W=sum_(i in S)w_i,
Q=(W,Q_sp), v=Q_sp/W and delta^2=1-|v|^2=Q^2/W^2.

For a temporal symmetric current H write A=H_sp and define

N(H,Q)=max( ||A||F, ||A*v||2/delta, |v^T*A*v|/delta^2 ).

At a singleton, delta=0 and transversality makes the last two numerators
zero; define the seed norm as the TT Frobenius norm. At generic larger
subsets delta>0.

The uniform bound for n=|S|:

N(H_S,Q_S) <= C_n * W^n/[product_i w_i] / kappa^((n-1)/2),

C_1=1,
C_n=64*sum_(set partitions pi of S, k=|pi|>=2)
      (k+1)!*32^(k+1)*13^k*product_(A in pi) C_|A|.

An explicit envelope is

C_n <=2*(2*10^10)^(n-1)*n!   for n>=1.

This controls every angular configuration and nested collinear approach
on the generic domain, with explicit soft-energy factors. It is NOT
energy-independent boundedness and does not assign values to excluded
exactly-collinear propagator points. Leaf norms multiply the RHS if
the leaves are not normalized.


## 5. Energy and multiplicity induction

Suppose the proposed bound holds on every child A. The energy product
in one root partition is bounded by

product_A [(W/W_A)*W_A^|A|/product_(i in A)w_i]
 =W^k*product_A W_A^(|A|-1)/product_i w_i
 <=W^n/product_i w_i.

Canonical powers combine to kappa^(-(n-1)/2). Multiplying the source
factor16 by the propagator factor4 proves exactly the C_n recurrence
stated above. The induction starts with unit physical TT leaves.

For c(z)=sum_(n>=1) C_n z^n/n!, labeled set-partition enumeration gives

c=z+2048*sum_(k>=2)(k+1)*(416c)^k
 =z+2048*((1-416c)^(-2)-1-832c).

Take cap=1/10^10 and radius=1/(2*10^10). Exact rational arithmetic shows
radius+2048*((1-416cap)^(-2)-1-832cap)<cap and416cap<1.
Iterating the nonnegative-coefficient fixed-point equation from c=z
stays bounded by cap at radius. Every coefficient stabilizes after
finitely many iterations because the nonlinear term begins at c^2.
Monotone convergence therefore bounds each coefficient by cap/radius^n,
giving C_n<=2*(2*10^10)^(n-1)*n!.

No infinite physical perturbative solution is inferred. This is a
numerical majorant for each finite tree coefficient; the majorant series
converges in its auxiliary counting variable.


## Scope of the auxiliary convergence argument

The convergent series counts a positive numerical majorant. It does not
construct an infinite interacting field expansion, exchange a physical
limit with a detector sum, or provide virtual corrections. Each asserted
physical source estimate is at a specified finite multiplicity.

The energy denominators are retained. A uniform angular bound with
explicit W^n/product(w_i) dependence is not a claim that the current
stays finite when any individual energy goes to zero with the others
fixed. Integrated soft subtraction and real-virtual completion remain
separate obligations.
