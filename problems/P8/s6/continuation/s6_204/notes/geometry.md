# Full four-dimensional geometric derivation

Write g=a^2 ghat with ghat=-deta^2+E_ij dx^i dx^j and E=exp(gamma). Because tr gamma=0, det E=1 exactly. The trace of the spatial extrinsic curvature of ghat is therefore zero, including its nonlinear terms. Its scalar curvature contains the time term tr[(E^-1E')^2]/4 and the intrinsic spatial scalar curvature.

At quadratic order, compact spatial integration gives

integral Rhat^(2)=integral[||gamma'||^2-||grad gamma||^2+2||div gamma||^2]/4,
Rhat^(1)=partial_i partial_j gamma_ij.

The scalar conformal shift is6a''/a. Its dependence on gamma vanishes exactly: det E=1, ghat00=-1 and the scale depends on eta alone. Hence sqrt(-g)R_old is a^2 times Rhat plus a gamma-independent term, while sqrt(-g)R_old^2=(Rhat+6a''/a)^2. Expanding these identities gives the full Einstein and R-squared quadratic actions stated in notes/conformal.md.

The Weyl background is zero and its quadratic density is the physical four-dimensional conformal Weyl quadratic form. The independent checks do not simply assume this simplification. They construct literal metric and inverse jets, all Christoffel symbols and derivatives, the full Riemann and Ricci tensors and scalar, and contract every index before taking the second metric-amplitude coefficient.

For a plane phase c=cos(P.x), s=sin(P.x), with independent tracefree matrix time jets M0,M1,M2, the exponential jets through second amplitude order include

E=I+epsilon c M0+epsilon^2 c^2 M0^2/2,
E'=epsilon c M1+epsilon^2 c^2(M0M1+M1M0)/2,
E''=epsilon c M2+epsilon^2 c^2[(M0M2+M2M0)/2+M1^2],
partial_i E=-epsilon s Pi M0-epsilon^2 cs Pi M0^2,
partial_eta partial_i E=-epsilon s Pi M1-epsilon^2 cs Pi(M0M1+M1M0),
partial_i partial_j E=-epsilon c PiPj M0+epsilon^2(s^2-c^2)PiPj M0^2.

All noncommuting products are retained. Exact inverse-product checks and the background/linear scalar-curvature coefficients are checked inside the independent fixture. Averaging the quadratic coefficients at phase quadratures averages c^2 and s^2 to1/2 and cs to0.

The resulting full Weyl quadratic coefficient has no conformal scale jets. For generic P and noncommuting Mj it is

||M2||^2/4-q||M1||^2+3||M1P||^2/2
-q tr(M0M2)/2+(M0P).(M2P)
+q^2||M0||^2/4-q||M0P||^2/2+(P.M0.P)^2/6.

Its relation to the compact Hessian includes a retained time divergence, detailed separately. This literal geometric calculation is a local-action check, not a substitution of a flat quantum state into the curved response.
