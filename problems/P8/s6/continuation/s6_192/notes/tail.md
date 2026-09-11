# Complete C2 current comparison and finite amplitude

Use the actual S6.191 weighted covariance error W=omega(Sigma-Sigma_ref).
Its full high-band integral norm bounds through parameter order2 are
T=(1e-65,1e-47,1e-28). Multiplying by the physical varying vertex gives

    integral |J_actual-J_ref| < 3 G0 T0 =6e-65,
    integral |partial_e(J_actual-J_ref)| <3(G1 T0+G0 T1)<1e-46,
    integral |partial_e^2(J_actual-J_ref)|
      <3(G2 T0+2G1 T1+G0 T2)<1e-27.

These retain the derivative of the observable, not just of its state.

The reference-minus-fourth-order integrand from contour.md has
radial measure bound (4/18)nu^2 dnu. Its complete infinite integrals
are below (4/18)B_a m^6 K^-2, hence below1e-13,1e-12,1e-10.
The finite-band fourth-order comparison has integral at most
B_a K^4/2. Adding this to the exact finite-band current estimates,
and adding both infinite tails, gives uniform bounds

    integral |partial_e^a(J_actual-J_ad4)| < (1e77,1e94,1e111)[a].

The partition is the same epsilon-independent proof band as before.
These are per-comoving-volume homogeneous currents for the actual
state. Both bands and all three constrained Proca modes remain.

At finite momentum the exact matrix flow and reference are smooth
in epsilon and t. The common integrable high-band bounds through
order2 and finite-band bounds justify differentiation under the
full integral, uniformly on the compact slab. For smooth detector D
the resulting current is C2 in epsilon into C0(I), up to the
endpoints of the admitted parameter interval.

The Banach-valued integral Taylor formula therefore proves

    ||Jcomp(epsilon)-Jcomp(0)-epsilon Jcomp'(0)||C0
        <= epsilon^2 * 1e111/2.

At epsilon0 both sides are zero; for nonzero epsilon the inequality
is strict because the uniform second derivative bound has a strict
margin. Dividing by the same fixed kappa=1e800 gives bounds
(1e-723,1e-706,1e-689), and Taylor coefficient1e-689. Dividing the
unimodular current by a^3>=1 to obtain proper density preserves
the bounds since a is epsilon-independent.

These small kappa-scaled numbers describe this mathematical
comparison only. Its finite matching to the original covariant
prescription is not yet proved. They are not a full physical
self-energy, inverse, finite-coupling solution or cutoff bound.
