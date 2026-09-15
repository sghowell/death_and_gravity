# Complete operators and homogeneous quantum forces

Fix either original cutoff c=1,2 and use the same physical seed, whitening
and reference symplectic/metaplectic evolution as before. Nothing below
assumes that the evolved state remains close to the vacuum or inside the
unmodified classical core.

## Definition with all scalar centers retained

In absolute homogeneous variables Y, let
h(u,Y,z)=kappa integral V Hbar(Nstar,Z) d^3x,
g=(h-h_ref)(u,Y,S_u z), g0=h(u,Y,0).
The fixed reference quadratic h_ref is independent of live Y and has
zero value at z=0. Set K_symbol=chi_c(g-g0). The complete interaction
Hamiltonians are
A_C=g0 I+Q_V[(1-D)K_symbol] and A_W=g0 I+OpW(K_symbol),
where D=Delta_w/4. The full scalar center is present in the energy and
the classical force. Factoring its phase later does not delete that force.

For F_Y=avg[exp(3v) R_full(u,Nstar)^(-3/4)], pulled back by S_u, set
F0=R_full(u,N0(Y))^(-3/4) and
F_ext=F0+chi_c(F_Y-F0).
Apply each of the two complete orderings separately. In particular F0
is not replaced by the clock value1 when differentiating live Y.
Physical volume is exp(3alpha) times this operator expectation.
Physical state and observable use the same fixed reference conjugation.

## Joint analytic bounds without an extra time assumption

The entire centered phase amplitude is bounded on the initial complex4R
phase ball and the homogeneous1e-122 polydisc. The full spatial estimate
of notes/geometry.md and the reference quadratic bound give
A_K<=2*1000*kappa*hvariation+1e128(8R)^2/2 <1e280.
The centered volume amplitude is<1e-510. Use these larger, fixed rational
ceilings A_H and A_F, not favorable sampled values.

The inherited simultaneous phase Cauchy radius R/8 and the complete
radial-cutoff derivative argument give, for0<=n<=196,
J_n(A)=2 A 2056^n (n!)^3/R^n.
Each ratio J_(n+1)/J_n is<1e-9. A homogeneous polydisc of radius
d=1e-124 fits strictly around every point of the real1e-390 path ball.
Cauchy differentiation in any one or two homogeneous coordinates
multiplies these bounds by j!/d^j, j=0,1,2 (the factor2 also covers
mixed second derivatives). These are complex PHASE/PARAMETER estimates;
the fixed profiles need only their original real-time C5 regularity.
The cutoff, reference flow and quantization are fixed in Y, so these
derivatives commute with them.

## Both operator orderings, not a symbol-norm shortcut

Use the independent normalized coherent-kernel Schur estimate proved
in S271, with its full constant1e112 and all phase coordinate orders<=2.
The exact heat-kernel integral remainder gives
||OpW(b)-Q_V[(1-D)b]|| <=1e112 *96^2 J_4(A)/32.
The coherent part obeys3A+96 J_2(A)/4. Thus
O(A)=3A+96 J_2(A)/4+1e112*96^2 J_4(A)/32
bounds BOTH full operators. The same expression with the homogeneous
Cauchy factors bounds every first and second parameter derivative.
All196 phase orders needed for the96-dimensional Schur estimate are
retained. A symbol supremum is not used as a Weyl operator norm.

The symbols are real at real Y and u. Their Weyl and coherent maps are
bounded self-adjoint operators, with norm-continuous time and parameter
dependence on this compact interval. Consequently they generate exact
unitaries; no Fock occupation projection is introduced.

## Complete canonical force budget

For k_i=<partial_i K> and D_h=kappa Vol a^3, all six corrections are
delta alpha'=k_p/(3D_h),
delta p'=(-k_alpha/3+pm k_pm+ph k_ph)/D_h,
delta pm'=-pm k_p/D_h,
delta eta'=1e100 k_ph/D_h,
delta ph'=(-1e100 k_eta-ph k_p)/D_h,
delta M1'=k_pm/D_h.
These are the absolute canonical chain rule proved in homogeneous.py,
not partial derivatives at a fixed density accidentally treated as
canonical derivatives. In particular a^3 pm remains conserved.

Since Vol>1, |a^-3|<2 and all density coefficients lie in the recorded
domain, the ceilings1e104 and1e106 dominate all five force-row sums and
their first parameter derivatives, including the heavy1e100 conversion.
With the complete operator bounds this gives
B_quant<=1e104 O(A_H/d)/kappa <1e-230,
L_quant<=1e106[O(A_H/d)+O(2A_H/d^2)]/kappa <1,
and5 max_i ||partial_i K|| <1e470.
Here O(A/d) denotes inserting the amplitude and the appropriate Cauchy
factor into the same linear bound, not a new allowed amplitude input.
These inequalities hold for every normalized state in the invariant
subspace. They do not assert a small norm for the state displacement.
