# Actual physical plane QEI to affine geometric curvature

All conventions and the physical reference come from the fully
replayed A.22/A.20 chain. Primes in this note are AFFINE lambda
derivatives. Since d lambda=a dt, a'=H and

    Hdot=a a'',
    Hddot=a a'a''+a^2 a''',
    Hthird=a a'^2 a''+a^2 a''^2+3a^2 a'a'''+a^3 a''''.    (3)

Choose the normalized compact spatial profile j_tau(z)=tau^-1/2
j(z/tau), where j is A.22's clamped polynomial. Its first derivative
norm is 3/tau^2 and integral j_tau j_tau'=0. The affine sampler v
is the cubic step chi from zero to one on [-tau/1000,0], followed
by h(lambda)=1-3(lambda/(2tau))^2+2(lambda/(2tau))^3 on [0,2tau].
All outer value and first derivative traces vanish, and the two
pieces have value one and derivative zero at their join.

Use f(t,z)=v(lambda(t))j_tau(z) in A.22, not a null-line sampler.
The physical plane measure is dvol=a dt dz=d lambda dz. Substituting
the actual proper-time derivatives into A.22's operators gives

    L/a=P j_tau,
    M=U j_tau',
    G/a=V j_tau+a^-2 v j_tau',
    P=a v''-2a'v'+2(a'^2/a-a'')v,
    U=a v'-2a'v,
    V=v'-2a'/a v.                                       (4)

After spatial integration the mixed quantum cross term vanishes
because integral j_tau j_tau'=0. The remaining quantum cost is

    Q=hbar/(8pi^2)[(5/9)integral P^2 d lambda
          +(1/3)(3/tau^2)integral a^-4 U^2 d lambda].     (5)

For an inhomogeneous state the analogous state cross term is NOT
discarded before applying the cap. Instead, w<=Phi_*^2 first bounds
the entire nonnegative G^2 weight, and only then can its purely
geometric cross term be integrated to zero. This gives

    W=Phi_*^2/3 [integral V^2 d lambda
                 +(3/tau^2)integral a^-4 v^2 d lambda]. (6)

Thus neither spatial homogeneity of w nor a lower bound on w is
used. The finite width tau is retained in both (5) and (6).

Applying (3) to the actual reference in A.22 yields

    Tconf_KK=hbar/(2880pi^2) N,
    N=-4a'^2 a''/a +beta_S[48a'^2 a''/a+84a''^2
                                       +72a'a'''+12a a'''']. (7)

Set r=-R_FK,KK=-2a''/a. The actual SEE gives r=kappa Ttotal_KK.
This TOTAL quantity is homogeneous because the geometry is FLRW,
regardless of homogeneity of the individual quantum state. On the
finite plane, Tother_KK>=ell and the scalar QEI give

    integral v^2 r d lambda
       >=-kappa Q-kappa W+kappa integral v^2 Tconf_KK d lambda
                              +kappa ell ||v||^2.        (8)

The spatial profile has norm one, so no area factor is missing.
Only one finite plane orientation is needed: FLRW makes the resulting
Ricci contraction the same on the outgoing rays used in the geometric
argument. This is NOT an equality of scalar stress on different
rays and does not impose state homogeneity or a global spatial cap.

The product sampler belongs to compact H0^2. Approximation by smooth
samplers from inside its affine-time and spatial rectangle gives
convergence in all required derivative norms after the smooth time
change. The relative stress and Wick expectations are smooth on
that compact region. The caps are imposed throughout this region,
so the same bounds persist in the limit. The geometry and field
are defined on a slightly larger open region when a putative affine
extension beyond 2tau is considered.
