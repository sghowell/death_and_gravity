# Entire state-difference response and full cutoff tails

All estimates use the full reference slab I=[-1/2,1/2]. The internal frequencies are nu=Omega_p and mu=Omega_(P-p). Directions are measured by v(D) from [hamiltonian.md](hamiltonian.md).

## Full pair products, not occupation-only terms

The exact identity delta V_D=(delta F_p)^t M_D F_q(T)+F_p(S)^t M_D delta F_q, together with the complete readout bounds, gives

|delta V_D|<=2*10^5 B sqrt(nu mu)(nu^-11+mu^-11)v(D).

It includes the error square through the actual T feature. Since the full S and T pairs are bounded by200² and400² times sqrt(nu mu)v, the ENTIRE two-time product difference obeys

|delta[V_D(t)conjugate(V_G(s))]|
<=10^11 B nu mu(nu^-11+mu^-11)v(D,t)v(G,s).

The reflected reverse term has the same bound. The i/2 factor accounts for both ordered branches; neither branch is omitted. The full volume product is at most64², so the complete memory integrand coefficient is below10^15 B.

For all external transfer, including large-small internal pairs,

mu<=nu+|P|/4<=L(P)nu,
nu<=L(P)mu,
L(P)=1+|P|/(4sqrt(n)).

With the full d³p/(2pi)³ measure,

integral Omega^-9=256/(105pi² n³)<n^-3,
integral Omega^-10=5/(8pi n^(7/2))<n^-7/2.

Both are exact full radial integrals, not cutoff integrals. Symmetry p<->P-p bounds the complete memory by10^16 B n^-3 L(P)v(D,t)v(G,s). The exact retarded triangle on a unit-length interval is bounded in L2_t by one; no external time derivative is needed and no upper endpoint is removed.

The entire one-mode contact difference is bounded by10^8 B n^-7/2 v(D)v(G), using ||delta F||( ||F_T||+||F_S||)<10^6 B Omega^-10 and a³/2<=32. The full fixed profile summand adds at most10^9 B n^-7/2 v(D)v(G), as proved in [profile.md](profile.md).

Their sum is below10^17 B n^-3 L(P). Define
norm_phys=||sqrt(1+|P|²)v(D)||L2(dt dP).
Since L(P)<=2sqrt(1+|P|²), the complete normalized matched response is bounded by

10^101/(kappa0 n³)

times the two physical norms. This is the whole state-selection difference plus its fixed profile, not the remaining comparison-state response.

## Full two-leg auxiliary cutoff

For K>=2sqrt(n), use the same reference projector Omega_p<=K on both memory legs and on the one contact leg. The fixed profile stays equal to its FULL unprojected value. This is an auxiliary mode projection, not a physical EFT cutoff.

Split external transfer into two regions. If |P|<=2K and either leg is missing, the reverse triangle inequality implies that BOTH internal frequencies exceed K/2. In the radial variable nu, p²dp<=64nu²dnu and pi²>9, so the entire ninth-power tail above K/2 is less than43K^-6. Its full memory bound is below10^17 B L(P)K^-6. Since K>=sqrt(n), K^-6<=n^-5/2/K.

If |P|>2K, retain the full untruncated integral bound and use
L(P)/(1+|P|²)<=1/K.
This consumes only the already declared spatial graph weights; it covers large-small pairs that a bare pointwise K^-6 assertion would miss.

The single-mode contact tail has radial Omega^-10 tail below K^-7 and is bounded separately. Combining the two external regions and this contact gives physical-graph normalized error below

10^102/(kappa0 n^(5/2)K).

## Actual common-clock graph and finite-cutoff contact

The first physical metric map gives Q=2(v+delta n_lapse)I and beta=iP B/a², delta=1/[2(1+t²)^3]. Thus v(physical)<=8sqrt(1+|P|²)(|n_lapse|+|v|+|B|). Both physical norms are controlled by eight times

norm_clock=||(1+|P|²)(|n_lapse|+|v|+|B|)||L2.

The full matched infinite-cutoff second-chart contacts cancel by the complete first-current identity. At FINITE K they do not cancel: the projected current's mean differs from the full fixed profile. The exact retained residual is given in [profile.md](profile.md). Its additional bound is10^9 B K^-7 and fits the declared allowance; it is not silently dropped.

The complete same-clock normalized bound is therefore

10^103/(kappa0 n³)<10^-1280,

and its full cutoff error is at most

10^104/(kappa0 n^(5/2)K)<10^-1180/K.

The actual n>10^196 and kappa0=10^800 prove these final strict displays. Both use the same stated graph and all external momentum. Neither estimate is a same-space inverse or nonlinear stability statement.
