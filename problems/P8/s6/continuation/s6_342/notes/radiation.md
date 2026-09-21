# Full selected physical TT amplitude

Use sum_i p_i+k=0, k^2=0 and J_i=epsilon(p_i,p_i)/(p_i.k).
All four external emissions evaluate the complete F at p_i+k, separately
for each i. For a triangle label, set P=p_a+p_b,v=P^2,u=(P+k)^2.

The complete amplitude is

 M5=1/(16pi^2 sqrt(kappa)) * {
   sum_i J_i [F(p_i+k)-F_sym]
   +g^2/4 sum_perm [A(v)*I_T(P,p_c,p_d)
       -2g^2 epsilon(P,P) Cbar(u;1,1)/((n-v)(n-u))]
   +g^4/4 sum_perm I_D(p_a,p_b,p_c,p_d) }.

I_T and I_D are the sums of all three and four scalar-line insertions.
Their explicit rational parameter densities are line_density in
radiation.py, integrated over all unit-cube coordinates. In a triangle
the heavy_angle coordinate is inactive, so its integration contributes1.
For every split line cyclically set Q0=0,Qj=-sum_(i<j)r_i.
With y=x0*gamma, the combined denominator is that of the N-line loop with
Q0 replaced by gamma*k. The exact-D stress insertion after the loop shift
is2epsilon(Q0-barQ,Q0-barQ), independent of y. The physical D4 normalized
density is

 -2 Gamma(N-1) x0 epsilon(Q0-barQ,Q0-barQ)
 /Delta_gamma^(N-1),

times the full simplex measure and S341 complex contour Jacobians.
Thus no internal-line weight or radial Gamma factor is omitted.

The common sign is fixed by the inverse variation and independently by the
known two-light bubble. Its two distinct line endpoints are
H*x/a*log(M_u/M_v) and H*(1-x)/a*log(M_u/M_v), a=P.k.
Together they give-2H*DD_B(v,u). The literal S295 heavy resolvent supplies
the same negative divided-difference sign. The outer triangle branch
multiplies Cbar(u), not Cbar(v): the loop receives P+k at that vertex.
For a test kernel C(v)=v, A(v)D_C+D_A*C(u)=D_(A*C);
the wrong endpoint fails at noncoincident original configurations.

## Leading soft subtraction and continuity

Let F0 be the Born hard kernel and S0 the Born soft current. Writing
F_i=F(p_i+k), the exact nonleading bracket is

 sum_i J_i(F_i-F0)
 +(sum_i J_i-S0)(F0-F_sym)
 +all_internal_and_outer_terms.

There is no replacement of the four off-shell hard coefficients by one
common value. The entire S341 complex gap holds along each split interval,
so internal terms are finite as omega->0. Equal invariants are treated by
the original gamma integral and its continuous divided-difference limit.
Only the external currents supply the leading1/omega pole.
No full unprojected Ward tensor or arbitrary curved completion is inferred.
