# Original local action, classical tensor tree and physical normalization

The homogeneous unit-Frobenius exponential tracefree amplitude is gamma. The original finite action is

Gamma_loc,2=(1/(64pi^2)) integral a^3[A gamma'^2-(D gamma)^2/60]dt,

A=5m^2/12-R0/36,
D=Dt^2+H Dt,
Ddag=Dt^2+5H Dt+2H'+6H^2.

The constant m^4 density has no unimodular tracefree variation. The first scalar-curvature variation vanishes, but the full R0-dependent second variation remains inA. These are the S189 coefficients of the actual fixed prescription, not a new subtraction.

Normalize the retarded current derivative byT_Gamma=64pi^2/a^3 times its gamma-coordinate derivative. The finite local part is

Tloc=-2[A L+A'Dt]-(1/30)Ddag D,
L=Dt^2+3H Dt.

The original canonical chart is gamma=2h/sqrt(kappa). There are TWO factors2/sqrt(kappa), one for the current and one for its variation. Consequently the canonical quantum force isT_Gamma h/(16pi^2 kappa). This agrees exactly with the S189 finite force: Tloc/(16pi^2 kappa)=-Qloc.

The actual classical tensor equation at zero transfer has Euler term-Lh. Thus the complete retained equation with canonical forcingf is

-Lh+T_Gamma h/(16pi^2 kappa)+f=0,
T_total h=-16pi^2 kappa f,
T_total=T_Gamma-16pi^2 kappa L.

Neither a canonical chain factor nor kappa has been removed. The zero-transfer rotational decomposition prevents linear coupling of this tracefree channel to the scalar clock/matter/constraint channels. No conclusion about their nonzero-transfer coupling follows.

After separating the original fourth local coefficient-1/30, the known finite-plus-classical local coefficients are

c3=-H/5,
c2=-2A-(4H'+11H^2)/30-16pi^2 kappa,
c1=-6HA-2A'-(H''+7HH'+6H^3)/30-48pi^2 kappa H.

There is no zeroth-order term from THIS explicit finite action and tensor tree. This does not set the complete current's remaining nonlocal or prescribed second-vertex contacts to zero. They are retained in the actualV_Gamma constructed in notes/remainder.md.

The code checks all five derivative coefficients by direct operator expansion, independently checks the finite-action Euler sign, and retains the complete variable-coefficient primitive kernels. The background, m1000, kappa1e800 and finite coefficients are unchanged.
