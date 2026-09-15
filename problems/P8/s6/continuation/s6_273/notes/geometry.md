# Complete spatial average and auxiliary estimates

All norms below are the inherited Wiener norms on the entire reconstructed
fields, not norms of a finite Taylor polynomial. The finite cutoff is on
input canonical phase coordinates only. Nonlinear harmonics and generated
means remain.

## Uniform joint domain and exact trace source

For |u|<=T=1e-180 and complex homogeneous distance<=1e-122 from Ybar,
a^2,a^-2,a^3,a^-3 have modulus<2. The background is spatially constant.
The fixed reference flow maps the initial4R phase ball inside8R; the old
field rows evaluated at8R are used without enlarging the input shell.
Recover v and tau from their scaled rows by4/(1+P)^2. The complete metric
and inverse A2 norms remain<2, so the inherited full determinant/shape
fixed point and formal-transpose mean-zero Neumann inverse apply.

A sharper bound than the old trace-momentum norm is essential here. For
pi_trace=(Pi_v/6) gamma^-1 and det gamma=a^6 exp(6v), direct contraction
with the FULL Lie derivative gives
pi_trace:L_xi gamma=Pi_v xi.grad(v)+Pi_v div(xi)/3.
Periodic integration by parts gives the exact trace generator
Dtrace=Pi_v grad(v)-grad(Pi_v)/3.
Only its tracefree projection vanishes; the fluctuation source does not.
geometry.py checks the Lie pairing for a completely nondiagonal metric and
all its first jets, not for a diagonal specialization.

With field rows f, p_b=8T+HOM and PV_b=100(T+HOM), use
Dtrace<=(PV_b+f_Piv) P v+P f_Piv/3,
Dshape<=180*32(1+P) f_Pitau,
Dvector<=12(1+P)^2 f_W f_PiW, and
Dmatter<=2[(1+f_deltaPiM) f_gradM+(2HOM+f_PiH) f_gradH].
The shape constant includes the entire DQ* and cotangent reconstruction,
and the vector term includes the full Gauss contact. The complete
adjoint inverse then gives lambda_A1<=4D0/P,
pi_TF<=5(32 f_Pitau+16 lambda_A1), shear<=192 pi_TF^2.
The spectral denominator P applies only to mean-zero modes; no inverse
zero momentum is taken.

## All twelve invariant images and their actual averages

geometry.bounds() records the complete pointwise images for
(p,G,delta_pm,ph,eta,shear,electric,magnetic,Wmass,M1gradient,Hgradient,
curvature) in their original source order. At input, every nonzero
coefficient has zero spatial mean, and the same free flow preserves its
wavevector. This does NOT imply zero mean after reconstruction.

Taylor's integral remainder for the entire density exponential gives the
recorded average bounds
|avg delta_p|<=10 p_b v^2+4 f_Piv v,
|avg G|<=36 v(1+P) f_PiW,
|avg delta_pm|<=10 v^2+12 f_deltaPiM v,
|avg delta_ph|<=10 HOM v^2+12 f_PiH v,
avg delta_eta=0.
All six following quadratic invariants retain their full pointwise
bounds in the average. The full curvature has linear term
-4a^-2 Delta v, which integrates to zero; its complete nonlinear remainder
is bounded by2e12 P^2(v+tau)^2. This includes nonlinear tensor curvature.
The inherited analytic metric/shape estimates give the same bound on
the complex phase domain; no sign cancellation is assumed.

## Whole lapse root and its integral remainder

All invariant images plus the homogeneous center remain strictly inside
1e-120. The source's500 full derivative bounds give
q<=1/30+1e9[rho_N+12(HOM+max image)]/3<1/20
on the NEW lapse ball rho_N=1e-115. At N=1 the full residual is bounded by
5e-400+12e4 HOM+1e4 sum(image).
The last term is essential: the self-map is for full phase-dependent
fields, not just the homogeneous center. The exact inequalities in
geometry.py give self-mapping and |C_N|>=3(1-q)>2.

The full implicit derivative formulas are
N_z=-C_z/C_N and
N_zw=-(C_zw+C_Nz N_w+C_Nw N_z+C_NN N_z N_w)/C_N.
The componentwise source envelopes imply |N_z|<=1e4 and |N_zw|<=1e18.
Every mixed contact is included. Taylor expansion is centered at the
actual homogeneous root N0(Y), not at an off-shell N=1. If Z and A are
the sums of the twelve pointwise and average bounds, respectively, then
|deltaN|<=1e4 Z and
|avg deltaN|<=1e4 A+1e18 Z^2/2.
These are exact integral-remainder bounds for the full implicit function.

## Full volume and reduced Hamiltonian averages

For U=R_full^-3/4, the same nonzero analytic branch gives
|U0|<2, |U_N|<12 and |U_NN|<1e4. Combining the exponential density and
lapse remainders bounds the entire centered average by
20 v^2+72 v |deltaN|+12 |avg deltaN|+5000 |deltaN|^2 <1e-510.
The formula retains the density/root cross term.

For the complete on-shell density f, the envelope identities are
f_z=H_z and f_zw=H_zw-C_z C_w/C_N. The respective bounds are1e3 and1e9,
with |f0|<100. Multiplying the entire exp(3v) gives the bound
10*100 v^2+1e3 A+1e9 Z^2/2
+6v(1e3 Z+1e9 Z^2/2) <1e-520
for avg[exp(3v)f]-f0. This does not replace the source by its quadratic
Taylor series; it bounds its exact remainder. The spatial measure and
homogeneous scale contribute at most2*1000*kappa to the Hamiltonian.
