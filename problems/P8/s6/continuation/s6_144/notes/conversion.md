# Literal parent conversion and the complete row bound

The MS amplitude minus the zero-outer-subtracted amplitude is
-F_MS sum_z C(z). The sign follows from the mixed Hessian variation;
as a check, varying the one-loop local bubble C^2/2 by a potential
quartic Gamma0 gives -C Gamma0. The zero-outer subtraction is local
in the parent with H retained: its reference counterterms are
delta L=3L F_D, delta g=g F_D, delta M=0 (delta G=G F_D/2).
These are assigned to this outer forest once.

At fixed common MS reference the matched star shifts are the negatives
of those finite values. Differentiating the actual tree amplitude
-L+g sum_z 1/(M-z) gives exactly -F_MS sum_z C(z).
Its b2 change is -2g F_MS/(M-2)^3, relative to the tree exactly -F_MS.
A local term in the parent is not necessarily a constant four-light
contact after eliminating H.

For each logarithm use the least integer n with argument<=2^n:
log(argument)<=n log 2<n. No floating logarithm is needed for the
certificate. Let n_m,n_M,n_T bound m^2,2M,T respectively, and put
V_c=L+g/(M-3), V=V_c+4g/M, P=24n_M^2+184n_M+184.
The four disjoint bounds are

 E_delta=2123366400000 N Y_hi^2 V/(Q_lo^2 sqrt(m)),
 E_bubble=3*3432 N Y_hi^2 V_c/Q_lo^2,
 E_triangles=3*64 N Y_hi^2 g P/(Q_lo^2 M),
 |F_MS| <=2 N Y_hi^2/Q_lo^2
             *[46+32n_m+(300+100n_T)/T].

The finite conversion contribution is tree_b2*|F_MS|.
At the actual reference the sum is approximately 1.84389805055489e-610,
or 4.60974512638723e-11 relative to the tree. The estimate is an
absolute bound and makes no amplitude-sign assertion.

This includes only this primitive's proper box forest and outer
reference conversion. First-order shifts in old loops, other finite
parameter/canonical-field changes, G1^2 and other order-two forests
remain separately owned; this calculation does not set them to zero.
