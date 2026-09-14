# Pinning the physical map to independent original coefficients

The frozen canonical parent defines, before solving auxiliaries,
U=R**(-3/4), Cchi=R**(-1/4), M=R**(1/4).
For a positive spatial conformal factor C with unchanged lapse, the
scalar temporal coefficient is C**(3/2)/N, its spatial coefficient
is N*C**(1/2), the Maxwell electric coefficient is C**(1/2)/N,
and its magnetic coefficient is N*C**(-1/2). Therefore scalar
and Maxwell coefficients independently require
C=Cchi**2=M**(-2)=R**(-1/2). The positive spatial density is
C**(3/2)=U, and every scalar mass or original scalar-source term
uses the spacetime density N*U.

These are independent coefficient families in the literal original
ADM action/Hamiltonian, not coefficients inferred from a volume
cancellation. At R=1 the actual C, U, Cchi and M all equal one.
The S261 expression R-1/2 instead gives C=1/2, volume sqrt(2)/4,
scalar gradient 1/sqrt(2), and Maxwell magnetic sqrt(2). It fails
even this reference metric and multiple independent action terms.

The full function R and all its derivative dependence remain present.
The declaration in source.py uses the immutable Cchi directly; the
separate M and U identities test that declaration. Neither a new
background nor a different action has been chosen to repair the error.
