# Independent graph polynomials and a uniform compact gap

For each refined graph let alpha_e>0 be its edge parameters.
Write Lalpha for their light sum and Halpha for their heavy sum.
Spanning-tree complements define

    U=sum_T product_(e not in T) alpha_e.

An independent weighted Laplacian cofactor, followed by
w_e=1/alpha_e and multiplication by product_e alpha_e,
reconstructs the same polynomial. Self-loops do not enter
the Laplacian but do remain in that prefactor. Every one
of the 192 comparisons is exact. U is positive in the
interior and homogeneous of degree two.

Spanning two-forest complements similarly define the degree-three
kinematic polynomials. A component containing one or three
external labels has invariant one (on-shell light mass squared).
Components with two labels define Ps, Pt or Pu. Empty or full
external-label components carry invariant zero. Momentum
conservation identifies each cut with its complement.

Set s=2+v, t=0, u=2-v, |v|<=1. The second polynomial is

    F=U(Lalpha+M Halpha)-Psingle-(2+v)Ps-(2-v)Pu.

There is no expansion of a heavy propagator. For M=16 define

    Fcoarse=U(Lalpha+16 Halpha)-Psingle-3(Ps+Pt+Pu).

The exact difference from F at that mass is
(1-v)Ps+3Pt+(1+v)Pu. Its real part is nonnegative throughout
the complex unit disc, because the cut polynomials have
nonnegative coefficients and -1<=Re(v)<=1.

The independent exact witnesses in notes/positive.md establish

    Fcoarse >= U(Lalpha+Halpha)/4.

The mass lift is exactly (M-16) U Halpha, hence

    Re F >= U[(Lalpha+Halpha)/4+(M-16)Halpha].

On the normalized simplex sum alpha=1, this is at least U/4>0.
The actual heavy mass squared is checked to exceed sixteen.
Thus none of these full bare graph denominators vanishes in the
strict simplex interior anywhere in the stated forward disc.

The derivative is dF/dv=Pu-Ps. Every pair-cut monomial occurs at
most once, since the complement fixes its spanning forest.
Connectedness supplies a crossing edge that extends that forest
to a spanning tree. Its monomial therefore occurs in U sum alpha.
Native coefficient checks confirm, for every refinement,

    Ps+Pt+Pu <= U sum alpha.

It follows that |dF/dv|/|F| <= 4 on the disc and interior.
This is a pointwise denominator-ratio bound only. It has not
been interchanged with a divergent bare loop integral.

On simplex boundaries U can vanish: the bare double bubble
at parameters (1,0,0,0) is an exact negative control. Neither
U/4 nor the derivative ratio proves convergence at those
boundaries. Ultraviolet subtractions remain necessary.
