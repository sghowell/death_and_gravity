# Explicit full-source Sobolev and Frechet bounds

All norms in this note are coordinate norms on the fixed CD slab,
not a change of physical metric. Let b=psi_t, c=psi_tt,
v_i=psi_i, w_i=psi_ti and D_ij=psi_ij. With p=du and a=(1+t^2)^2,

    X=(1+b)^2-|v|^2/a^2,
    Boxu=c+3H(t)(1+b)-tr(D)/a^2,
    Z=(1+b)^2 c-2(1+b)v.w/a^2
           +v.D.v/a^4+H(t)(1+b)|v|^2/a^2.

The last expression is also computed from the full covariant Hessian;
its connection terms must be retained. Since a>=1, |H|<=2, |H'|<=4
and every psi jet through3 has magnitude at most delta<=1/100,

    |X-1|, |partial_i X| <=3delta,  ||p||Euclidean<=2,
    |E|<=40delta, |partial_i E|<=42delta,
    |Z|, |partial_i Z|<=2delta.

For example X-1 is bounded by2delta+4delta^2, and its spatial
derivative by2delta+8delta^2. For E=Boxu-3H(u)X,
write Boxu=3H(t)+Box psi and split H(t)-H(u), bounded by4|psi|.
This gives4+6+12+18=40 for the E coefficient. Differentiation gives
10+12(11/10)+18<42. For Z and its spatial derivative, the exact
displayed polynomial gives respectively

    [(1+delta)^2+12(1+delta)delta+9delta^2]delta <2delta,
    [(1+delta)^2+26(1+delta)delta+39delta^2]delta <2delta.

Thus X stays between97/100 and103/100, strictly within the coefficient
strip. The complete R bounds give
|R-1|<=3delta and |partial_i R|<=4delta.

For Q=E/X+3R_u/(4R)+C Z, source.md gives |C|<6,
|C_u|<15 and |C_X|<32. Hence |Q|<=70delta. The spatial derivative
has three contributions. Their coefficients, divided by delta, obey

    E/X:       42/(9/10)+120delta/(9/10)^2,
    3R_u/4R:  81/5+(351/2)delta,
    C Z:       12+222delta.

Their sum is below82 at delta=1/100. These bounds use the full
R_uu,R_Xu,R_XX, not derivatives of a truncated source.

To obtain L2 bounds WITHOUT a support-volume estimate, let
G_psi^2=sum_(|alpha|<=3)|partial^alpha psi|^2. Every individual jet
is at most both delta and G_psi. The linear defect estimates above
therefore also hold with the final delta replaced by G_psi.
Apply the product rule to S=p(R-1)Q:

    ||S||pointwise <=420delta G_psi <=512delta G_psi,
    ||partial_i S||pointwise
       <=[420delta+560+492]delta G_psi
       <=2048delta G_psi.

The full source value also has the bound420delta^2. Spatial integration
gives ||S||2<=512delta U_psi and
||partial_i S||2<=2048delta U_psi, U_psi=||G_psi||2.

## An arbitrary variation, not only a spatial derivative

For a compact smooth test eta, vary u at FIXED metric. The exact
derivative of X is2u^mu partial_mu eta, bounded by3G_eta.
The corresponding bounds are

    |D R[eta]|<=4G_eta, |D E[eta]|<=42G_eta,
    |D Z[eta]|<=2G_eta, |D Q[eta]|<=82G_eta.

For D Z, differentiate the displayed exact polynomial in b,c,v,w,D;
the same coefficient1+28delta+66delta^2 is below2.
For D E include -3H'(u)eta X-3H(u)D X and the full Box eta.
The Q derivatives use the same complete coefficient bounds as above.
The three terms in D[p(R-1)Q] then give

    ||D S[psi]eta||pointwise <=2048delta G_eta,
    ||D S[psi]eta||L2 <=2048delta U_eta.

Only the background psi jets need be small; eta is an arbitrary
compact smooth direction. H3 contains all coordinate derivatives
through3. It is a conservative test norm, not a claim that every
direction or metric perturbation satisfies the same estimate.
