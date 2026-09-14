# Uniform nonlinear root, not a linear implicit-function assumption

On |N-1|<=10^-6 and every stated invariant coordinate |z_i|<=10^-250,
the whole independent original source-jet boxes are evaluated with
exact outward rational intervals. Their R bounds stay inside
[999/1000,1001/1000]. The literal constraint obeys

    -31/10 < C_N < -3,
    |C_NN|<1000, |C_NNN|<10000.

The reported rational outward endpoints are rounded outward from
the exact intervals, never converted through binary floats.
Conservative integer ceilings give
sum_i|C_z_i|<100 and an independent full source-jet gradient sum<1000.

The center residual is evaluated by a CENTERED mean-value argument.
The bare comparison C(1,0)=0 is exact. Move all independent source
atoms from their bare center values to their full values along a
segment contained in the verified positive-R jet box, and then
move the twelve parameters from0 to z. If g_i and l_j are the
reported outward gradient ceilings, this gives

    |C_full(1,z)| <= b
       := sum_i g_i*10^-250+sum_j l_j*10^-380 <10^-245.

The smaller source-only bound applies at z=0. Every entire
normalized-source error is below the common10^-380 bound, so none
is dropped. The heavy mass ratio is held at its actual value in the
proved interval; the bare zero identity is independent of it.
This avoids subtracting separately rounded O(1) fourth roots, which
would inherit the older interval engine's100-bit enclosure floor.
The frozen engine is unchanged.

Let F_z(N)=N+C_full(N,z)/3. On the outer strip
|F_z'|<=1/30. For rho=10^-245,

    |F_z(N)-1| <= b/3+rho/30 <rho.

Thus the complete closed inner interval maps strictly into itself.
The contraction theorem gives an actual full-source root, not a
finite Taylor-series solution. Strict negativity of C_N proves
uniqueness on the whole outer strip. No uniqueness outside it is
claimed. With N0=1 the error after k iterations is bounded by
10*b/(29*30^k). The mean-value theorem also gives |Nstar-1|<=b/3.
At z=0 only the source-only part remains, giving the stated
10^-375 bound, without deciding whether the shift is zero.

All input variables lie in a finite-dimensional analytic extension
of the physical invariant data. Restriction to physical metric,
matter and vector jets preserves this local root conclusion.
It does not show that an arbitrary invariant tuple satisfies the
spatial momentum constraint, gives Cauchy data for nonlinear
evolution, or is the outcome of a joint quantum measurement.
