# Symbolic-dimensional geometry and contractions

The continuation is the literal one in FORMULATION.md. It keeps the
four-dimensional invariant coefficients and determinant exponent,
while every index sum and determinant uses D+1 dimensions. It is a
new choice only for the previously unspecified Y^2 derivative terms.

## Curvature from a warped fiber

Use a two-dimensional base metric diag(N^2,A^2) in (u,x), with a
flat (D-1)-dimensional fiber of scale A. Put n=log N,
v=log(A/a), L=H+v_u. The independent Ricci components are

    R00=-D(L_u+L^2-n_u L)
        -(N^2/A^2)[n_xx+n_x^2+(D-2)n_x v_x],
    R01=(D-1)(n_x L-v_ux),
    R11=-(A^2/N^2)[L_u+L(D L-n_u)]
        -n_xx-n_x^2+n_x v_x-(D-1)v_xx,
    Rii=-(A^2/N^2)[L_u+L(D L-n_u)]
        -v_xx-(D-2)v_x^2-n_x v_x, i>=2.

For the Riemann square use the base sectional curvature K,
the base Hessian C_ab=(nabla_a nabla_b A)/A and
S=|grad A|^2/A^2:

    Riem^2=4K^2+4(D-1)|C|^2+2(D-1)(D-2)S^2.

All components and contractions at D=3 agree with S6.61's full
coordinate Christoffel/Riemann calculation. The formula is a
polynomial identity in D, not a fit to integer-dimensional samples.

For temporal and spatial deviations t,s, the continued auxiliary
metric and relative measure are

    N^2=(1+s)^(D/2)(1+t)^(-1/2),
    (A/a)^2=(1+t)^(1/2)(1+s)^((D-2)/2),
    sqrt(g_tilde)/sqrt(g)=(1+t)^((D-1)/4)
                         (1+s)^(D(D-1)/4).

These powers, the derivative jets and the full scalar heat
coefficient are expanded algebraically to quadratic mass order.
The opposite-momentum polarization is performed directly on
quadratic monomials, avoiding unnecessary higher-degree expansions.

## Exact Einstein sums

There are two distinguished frame indices: time 0 and the spatial
momentum axis 1. Each dummy label is assigned to one of them, an
existing anonymous transverse class or a new anonymous class.
A pattern using r distinct transverse classes has multiplicity

    (D-1)(D-2)...(D-r).

Summing these canonical patterns counts every Einstein assignment
exactly once. The total assignment weight for k labels is
(D+1)^k. The executable checks verify k=1,2,3,4, every individual
two-derivative invariant and their full sum at D=3.

For a mass tensor differentiated at most twice, a component has
rank at most four. Reflection symmetry requires every anonymous
transverse index to occur an even number of times. At most two
distinct anonymous indices remain, so their value is exactly
represented by the frozen four-dimensional frame after a local
relabeling. Dimension dependence resides in multiplicities and
Ricci/scalar contractions, not an unproved interpolation.

The continued operators are homogeneous of degree two in Y in
every D. Their value and first variation at Y=0 vanish identically.
The prior flat potential, linear counterterms and selected stress
profiles are not changed.
