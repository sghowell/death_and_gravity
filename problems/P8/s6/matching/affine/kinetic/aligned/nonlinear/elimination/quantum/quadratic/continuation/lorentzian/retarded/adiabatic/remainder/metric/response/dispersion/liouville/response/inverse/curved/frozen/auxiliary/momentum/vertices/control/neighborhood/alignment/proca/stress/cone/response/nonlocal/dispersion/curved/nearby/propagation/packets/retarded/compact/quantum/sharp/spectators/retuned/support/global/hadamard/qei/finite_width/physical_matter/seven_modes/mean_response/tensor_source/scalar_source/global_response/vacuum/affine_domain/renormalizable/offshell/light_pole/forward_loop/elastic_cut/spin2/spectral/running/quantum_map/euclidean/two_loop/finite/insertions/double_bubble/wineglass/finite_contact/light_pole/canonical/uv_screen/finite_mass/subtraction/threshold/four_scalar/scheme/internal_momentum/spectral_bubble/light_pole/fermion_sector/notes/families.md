# Noncommuting scalar expansion and the nine primitive rows

After exact H elimination write the scalar Hessian
S_B''=D^-1+W2, where W2 is homogeneous degree two in the
Phi background. It is the full nonlocal Hessian of the
actual quartic action, not the local L term alone.

Let F0,F2,F4 be the background-degree coefficients of
the scalar Hessian of the entire fermion determinant.
The subscripts denote background degree, not derivative
factorials to be supplied later. Formally,

    (D^-1+W2)^-1=D-D W2 D+D W2 D W2 D+...
    F''=F0+F2+F4+...

Keeping the operator order gives scalar terms

    degree 0: (1/2) Tr[D F0],
    degree 2: (1/2) Tr[D F2-D W2 D F0],
    degree 4: (1/2) Tr[D F4-D W2 D F2+D W2 D W2 D F0].

The two-point terms have coupling classes Y^2 and
Y times one tree quartic/heavy kernel. The four-point
terms have classes Y^3, Y^2 times one such kernel, and
Y times two such kernels. These are operator families,
not a claim that each is a single ordinary Feynman graph.

At A=0 the tree gauge covariance is Phi-independent.
Its three contributions through degree four are
(1/2)Tr[D_A F_AA,0], (1/2)Tr[D_A F_AA,2] and
(1/2)Tr[D_A F_AA,4], with coupling classes a, aY and
aY^2. Together the scalar and gauge expansions give
two vacuum rows, three two-point rows and four
four-point rows.

Setting P=-D F0 D, the variation of the scalar
quadratic half-trace is -(1/2)Tr[D W2 D F0],
the paired sector of S6.136. Varying
-(1/4)Tr[D W2 D W2] gives
+(1/2)Tr[D W2 D W2 D F0]. S6.135 bounds only this
last family's local-quartic vertex subset. The other
heavy-kernel pieces do not inherit that bound.

The expansion is formal through four derivatives at
the vacuum, with all momentum dependence retained.
It is not a finite-field convergence proof, a
hard-region-only matching calculation, or a complete
two-loop error budget.
