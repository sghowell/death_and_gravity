# Actual timelike-plane positivity and physical improvement

The plane embedding i(t,z)=(t,0,0,z) has purely transverse spacelike
conormal covectors. These cannot meet the nonzero null covectors in
a Hadamard two-point wavefront set. Therefore i x i pulls back the
two-point distribution, its real differential operators and its
positive type. Positivity can equivalently be checked by approximating
the transverse delta distributions with smooth test functions before
taking the well-defined distributional restriction. The target minus
vacuum two-point function is smooth and symmetric and has a smooth
diagonal after restriction. This argument is unavailable for the
vacuum two-point distribution pulled back to a single null line.

No transverse field average is squared: we restrict the expectation
of the renormalized local four-dimensional physical stress itself.
For any compact real f, the full improvement gives the identity

    f^2 T_DD = (1-2xi) f^2(D Phi)^2
              +2xi [D(f Phi)]^2-2xi (Df)^2 Phi^2
              -2xi D(f^2 Phi D Phi).                     (5)

The last term integrates to zero on the plane. The same identity
holds for smooth normal-ordered two-point differences, so no product
of unsmeared quantum fields or divergent boundary term is taken.
The first two terms have nonnegative coefficients on 0<=xi<=1/2.

To bound them, form the two compact real-operator kernels

    K1(x,x')=f(x)f(x') D_x D_x' W(x,x'),
    K2(x,x')=D_x D_x'[f(x)f(x')W(x,x')].                  (6)

Each is of positive type for a positive W. Their target-reference
differences are smooth, symmetric and compactly supported. Fourier
inversion on the diagonal identifies the integral of each difference
with the integral of its Fourier diagonal over all alpha, with factor
1/(2pi)^2. Symmetry under alpha->-alpha makes this twice its integral
over alpha_t>=0, with ALL alpha_z retained. Dropping the nonnegative
target term gives minus the reference Fourier-diagonal integral.
Smooth cutoffs on this half-plane justify the inequality before the
limit; the explicitly finite reference integral below justifies the
limit. The compact smooth difference has rapid Fourier decay.

Our Fourier transform of f uses exp(i theta_t t+i theta_z z). The
actual four-dimensional vacuum, restricted to the plane, is

    W0= hbar/[2(2pi)^3] integral d^3k/omega
                            exp(-iomega Delta_t+ikz Delta_z)
       =hbar/(8pi^2) integral_{omega>=|kz|} d omega d kz
                            exp(-iomega Delta_t+ikz Delta_z). (7)

The angular transverse measure is 2pi p dp and p dp=omega d omega.
Thus the mass-shell factor is not that of a two-dimensional field.
In (6)'s reference diagonal, put q=omega-kz and
theta_t=alpha_t+omega, theta_z=alpha_z-kz. K1 contributes
q^2 |fhat(theta)|^2. Integration by parts in K2 contributes
(alpha_t+alpha_z)^2 |fhat(theta)|^2, not q^2 or (theta_t+theta_z)^2.
Let T=theta_t, Z=theta_z; then alpha_t+alpha_z=T+Z-q.
At fixed theta the remaining region is 0<=omega<=T, |kz|<=omega,
empty unless T>=0. Its exact integral is

    integral_0^T d omega integral_-omega^omega d kz
            [(1-2xi)q^2+2xi(T+Z-q)^2]
       =(2/3)(1-xi)T^4+(4xi/3)T^3 Z+2xi T^2 Z^2.         (8)

This nonnegative expression is integrable against every smooth
compact sampler's rapidly decreasing Fourier transform. The
reference prefactor is hbar/(8pi^2) times 2/(2pi)^2. Because f is
real and the polynomial (8) is invariant under (T,Z)->(-T,-Z),
half-plane Parseval contributes (2pi)^2/2. In particular the T^3 Z
term does NOT vanish for a general sampler. It becomes <f_tt,f_tz>.
Equations (5)-(8) prove (1) with its stated normalization.

Writing L=f_tt and M=f_tz, the derivative density equals

    2xi(M+L/3)^2+(2/9)(3-4xi)L^2,                       (9)

nonnegative on the entire stated coupling interval. At xi=0 its
cost is hbar||f_tt||^2/(12pi^2), the minimal known-answer control.
The remaining state weight 2xi(Df)^2 is nonnegative pointwise.
Hence a ONE-SIDED cap w<=Phi_*^2 suffices to give the displayed
state-cap bound; an absolute cap and a physical momentum cutoff
are not inferred. All mass-shell, polynomial and prefactor
identities are checked by native algebra in flat.py.
