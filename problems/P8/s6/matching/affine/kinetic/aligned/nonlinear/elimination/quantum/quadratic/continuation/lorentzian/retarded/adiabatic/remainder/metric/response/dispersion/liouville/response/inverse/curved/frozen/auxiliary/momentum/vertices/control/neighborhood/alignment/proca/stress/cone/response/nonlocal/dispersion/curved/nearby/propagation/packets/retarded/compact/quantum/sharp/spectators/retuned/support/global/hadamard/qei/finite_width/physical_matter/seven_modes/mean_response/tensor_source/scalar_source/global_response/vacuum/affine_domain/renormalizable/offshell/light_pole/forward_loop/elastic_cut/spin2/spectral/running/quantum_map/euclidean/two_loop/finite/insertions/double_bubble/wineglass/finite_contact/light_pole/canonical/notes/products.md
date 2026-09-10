# Disjoint products and the fundamental cubic vertex

Choose all vertex-disjoint pairs of logarithmic proper cores in the
4608 raw graph records and contract both simultaneously. There are
72 such pairs and 12 distinct labelled tree cographs. Overlapping
proper cores are not admitted as pairs.

Independently start from the full heavy-exchange tree: two cubic
vertices, one heavy line, six external-label assignments, weight 1/2
per assignment. Construct three possibilities:

- replace both cubic endpoints by first-order cubic counterterms;
- replace either one endpoint and insert one heavy mass counterterm;
- insert two identical heavy mass counterterms in series.

The last case has the ordinary Dyson coefficient one, not an extra
factor two. All 12 graph keys, absolute Wick weights and signed
coefficients match the contracted two-loop forests. The three channel
sums are individually checked, preserving s/t/u assignments.

Let G1=LG Ibar/2 and M1=G^2 Ibar/2. The complete product
amplitude per channel h=1/(M-z) is

    G1^2 h - 2G G1 M1 h^2 + G^2 M1^2 h^3
      = Ibar^2 [L^2 g h/4 - L g^2 h^2/2 + g^3 h^3/4].

The h term is as important as the iterated heavy-mass h^3 term.
It comes from two endpoint cubic counterterms. It is not a local
momentum-independent quartic contact. Its diagnostic forward
second coefficient is 2G1^2/(M-2)^3, not zero.

For fundamental canonical vertex parameters
L_total=L+ell L1+ell^2 L2,
G_total=G+ell G1+ell^2 G2,
M_total=M+ell M1+ell^2 M2,
literal expansion of -L_total+G_total^2 sum_z(M_total-z)^-1
gives the order-two tree term

    -L2 +(2G G2+G1^2) sum h
        -(G^2 M2+2G G1 M1) sum h^2
        +G^2 M1^2 sum h^3.

The inherited new-order family quantities labelled delta g are
linear local tree variations 2G G2 for those families. Their
disjoint first-order products are already represented by the
separate forest terms above. This is consistent with their stated
family-only scope and does not revise an ancestor.

If instead one parametrizes the total squared coupling directly,
its second-order coefficient is delta g_total,2=2G G2+G1^2.
One must not add the cubic square again after using that total.
The native symbolic expansion checks both parametrizations.
All G1/G2 here are total CANONICAL vertex counterterms, not
bare coupling coefficients before canonical field factors.
