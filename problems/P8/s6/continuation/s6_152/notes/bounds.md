# Full-propagator bounds and the exact quadratic coefficient

On the inherited radius-one forward disc, |z|<=3,
|C|<2L, |I_R|<2/Q and

    |T_A+T_B|<16g log(4M)/(QM).

These are all-radius bounds on complete heavy propagators.
Using g/M<L/3, the three-channel linear-scale difference has

    E_linear < [24+16ellH] ell L^3/Q^2,
    ellH=log(4M).

The factors are respectively
(3/4)*2*(2L)^3*(2/Q)*ell/Q and
(3/4)*(2L)^2*16g ellH/(QM)*ell/Q.
Cauchy's estimate at radius one bounds its second forward
coefficient with no extra factor 1/2.

For the quadratic term, use its exact forward coefficient
instead of bounding a large constant amplitude. Put D=M-2,
s=2+nu, u=2-nu and t=0. Direct differentiation gives

    b2[sum_channels C^3]/b2_tree
       =3(-L+g/D)(-L+2g/D).

The actual 0<g/D<L/2 implies the absolute value is <=3L^2.
Therefore E_square/b2_tree<=3ell^2 L^2/(4Q^2).
The constant t-channel is retained but has zero nu derivative.

With ell<1200, ellH<600 and Q>144, the linear and quadratic
bounds are approximately 6.67799149312032e-612 and
1.09139364212751e-1007. The old-family bound is
19224 L^3/20736, about 1.11161080340594e-614.
Their exact rational sum obeys

    E_double_bubble,MS<6.68910760115439e-612<7e-612,
    E_double_bubble,MS/(4lambda)<1.67227690028860e-12<2e-12.

The decimals are diagnostic; rational comparisons prove the
bounds. No physical loop-momentum cutoff or heavy-mass series
was used, and this fixed-order estimate is not a higher-loop
truncation theorem.
