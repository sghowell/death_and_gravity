# Explicit finite Weyl ordering theorem and boundary

## Unchanged object

Retain the entire S6.270 source and its report SHA-256
5b6adbb4be256ec2df0051d781d43a2e94b8ab12d099c79295e87c821fb0fc01.
There are d = 48 unit-CCR pairs and m = 96 real phase coordinates after the
actual symplectic whitening of the original pure covariance. Keep R = 10^20,
T = 10^-2000, the full free reference, scalar phases, original matter/vector
state, residual translation constraints and the complete nonlinear spatial
and auxiliary reconstruction.

At every real |u| <= T, the full interaction g and normalized physical volume
F, already pulled through the reference flow, are holomorphic on the complex
initial phase ball 4R, with |g| < H = 10^1000 and |F - 1| < E = 10^-255.
Only phase analyticity is used here; the original real-time profiles are not
silently upgraded from C5 to analytic or infinitely differentiable profiles.

Let b(y) = exp(-1/y) for y > 0 and zero otherwise, and
theta(y) = b(1-y)/(b(y)+b(1-y)). Keep precisely

    chi_c(w) = theta((|w|^2 - R^2)/(c R^2)),  c in {1, 2}.

Both cutoffs equal one on the core ball R and vanish before radius 2R.
For f = g or F, set f_ext,c = f0 + chi_c(f - f0), with f0 = f(u,0).
The actual time-dependent scalar f0 is never reset. Write D = Delta_w/4.

## Operators being compared

The established calibrated-coherent operators are

    A_C,c = Q[(1-D)g_ext,c],   F_C,c = Q[(1-D)F_ext,c].

This continuation defines the distinct finite Weyl ordering

    A_W,c = OpW(g_ext,c),      F_W,c = OpW(F_ext,c).

The same unit-CCR Weyl convention, seed and symplectic whitening are used.
The identity Q(a) = OpW(exp(D)a) is the previously normalized same-seed
coherent identity, not an assumption that fixed-window coherent quantization
is covariant under arbitrary symplectic changes.

## Established estimates

A direct Gaussian-frame product-Schur proof gives, at this fixed dimension,

    ||OpW(a)|| <= 10^112 max_{alpha_j <= 2} ||partial^alpha a||_infinity.

All indicated mixed derivatives, through total order 192, are included.
The full cutoff/Cauchy amplitude jets through order 196 obey

    J_n(A) = 2 A 2056^n (n!)^3 / R^n,  0 <= n <= 196,

for A = H or E. Their successive ratios are less than 10^-9.
For the heat remainder r_f = exp(D)(1-D)f_ext - f_ext this implies

    ||OpW(r_f)|| <= 10^112 (96^2/32) J_4(A).

Consequently the volume ordering error is less than 10^-200, and the
Hamiltonian ordering error is less than 10^1056. Since
||F_C,c - I|| < 2E, the concrete Weyl volume satisfies

    F_W,c > I/2,          ||F_W,c - I|| < 10^-199.

Both complete real bounded Hamiltonians generate exact unitaries. Unitary
Duhamel comparison, with the entire scalar phases retained, gives

    ||U_W,c(u) - U_C,c(u)|| < 10^-944,
    ||U_W,c(u) - I|| < 10^-943.

On the original normalized seed, changing both state and readout between the
two orderings changes the volume mean by less than 2 * 10^-200.
Comparing c = 1 and c = 2 through the established coherent comparison gives
Weyl-state difference less than 10^-943 and Weyl-volume-mean difference less
than 10^-199. The original positive coherent POVM assigns the evolved Weyl
state outside-core probability less than 10^-942. This is not exact support.

The full translation symmetry and common zero-charge reducing sector remain.
Physical-picture comparisons conjugate both states and corresponding readout
operators by the same original reference propagator.

## Non-claims

This is neither generic positivity-preserving Weyl quantization nor an
original interacting volume mean. It does not construct the singular
unregularized Hamiltonian, remove either regulator, match physical Wilson
coefficients, bound omitted loops, establish all-energy UV positivity or
prove nonlinear global completeness. None of the constants is asserted
uniform in the mode count, torus size, cutoff choice or a longer time interval.

The stronger S6.270 coherent-only cutoff estimates are not reassigned
unchanged to Weyl ordering. Original V/G/B/P8 remain OPEN. The certificate is
an exact source-bound computational audit of written mathematics, not a
formal proof assistant certificate for its analytic or operator arguments.
