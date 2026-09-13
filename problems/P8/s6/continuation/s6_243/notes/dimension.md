# Full six-invariant dimensional geometry

Keep the six physical bilinear invariants fixed:

T=tr(DG), V=e^t DG e, W=(eDe)(eGe),
S=(tr D)(tr G), UD=(tr D)eGe, UG=(tr G)eDe.

For a unit internal direction nvec=u e+sqrt(1-u²)vvec, vvec lies on the transverse unit sphere in b=d-1 dimensions. Its second moment is delta_ij/b and its fourth moment is (delta_ij delta_kl+delta_ik delta_jl+delta_il delta_jk)/[b(b+2)]. Odd transverse moments vanish. Contracting these tensors with the FULL symmetric D and G gives, with R=1-u²,c=u(y-u),y=P/r,

geo00=S,
geo01=-R S/b+(R/b+c)UD,
geo10=-R S/b+(R/b+c)UG,

geo11=R²(2T-4V+S-UD-UG+3W)/[b(b+2)]
+R(y-2u)²(V-W)/b
-Rc(UD+UG-2W)/b+c²W.

This is a direct arbitrary-d tensor derivation, not interpolation from integer dimensions. Distinct trace/gradient orders are kept through every source time derivative. The remaining sphere average is
average u^(2j)=(1/2)_j/(d/2)_j, average u^(2j+1)=0.
Multiply the y-polynomial by the endpoint x-polynomial BEFORE this average. This retains all four ordered pair products and all mixed powers.

The expressions are analytic near d=3: the transverse denominators have no pole on |d-3|<=1/4, and only finitely many displayed moments occur. At complex d, sharp is algebraic and keeps d fixed. No positive noninteger-dimensional Hilbert norm is asserted.

Independent literal Cartesian moments in dimensions3,4,5,6 test two noncommuting matrix pairs, all four ordered geometries and every extra power u^0..u^4:160 exact comparisons. The direct Cartesian rule integrates every monomial using products of odd double factorials divided by d(d+2)...; it does not invoke the invariant formula.

For the covariant check use S216's generic full exponential-metric Hessians H_R,H_R2,H_Ric2,H_Riem2, including a^(d-2),a^(d-4), all trace/volume terms and ordered time integrations. Only the GEOMETRIC identities are reused. The scalar pole is derived with the S240 weights,

H_pole=-m² H_R/3+H_R2/36+(H_Riem2-H_Ric2)/90.

The m^4 volume Hessian is independent of P and cancels in this spatial difference, not in the whole determinant. All six invariant coefficients satisfy H_pole|d3=16 sum_r f4r|d3 Gamma_r, where f4r is the complete logarithmic radial coefficient. The old vector weights do not satisfy this scalar identity.

This proves the local continued coefficient and pole input. A dominated limit of the complete subtracted momentum/time integral requires an additional envelope for the full kernel; analyticity of these finitely many coefficients does not supply it.
