# Uniform actual W8 symbol limit at fixed external momentum

Let x = 1/r and Omega = omega/r = sqrt(a^-2+m^2 x^2). The complete eighth-order frequency has the scaled expression

W8/r = Omega + sum_(j=1)^4 P_j(t,z) x^(2j) Omega^(1-2j),
z = 1/(1+m^2 a^2 x^2).

The actual P_j and their required total time derivatives are bounded on the fixed compact CD slab near z = 1. The scale factor is positive there. Consequently the scaled expressions, their required time derivatives, and the normalized field readouts are smooth at x = 0, with uniform compact-time bounds. This uses all four W8 coefficients, not a frozen-frequency replacement of the full endpoint.

The scaled total time derivative includes the derivative of Omega and all total derivatives of P_j. It tends to -H/a. Thus the H and W'/W terms in the normalized momentum are multiplied by x and vanish in its leading limit. Exact checks retain these terms before taking x = 0.

At fixed external P, l/r = -n+P/r, and |l|/r = 1+O(1/r) uniformly in n. Real polarization frames can be chosen in a finite angular atlas. Their complete transverse projectors and longitudinal outer products are frame independent, so no global smooth individual frame is assumed.

The actual summed endpoint therefore has E0 = r h + O(1), uniformly in angle and on the fixed time slab for fixed P and tensors. On the lost shell r is between K-|P| and K. The O(1) error integrates to O(K^2) at fixed P, which is o(K^3). The exact leading shell calculation then gives the coefficient stated in FORMULATION.md.

This is a fixed-P symbol asymptotic. No bound uniform in all external momenta, no new Sobolev response estimate, and no analyticity of the prepared Borel-state data are claimed.
