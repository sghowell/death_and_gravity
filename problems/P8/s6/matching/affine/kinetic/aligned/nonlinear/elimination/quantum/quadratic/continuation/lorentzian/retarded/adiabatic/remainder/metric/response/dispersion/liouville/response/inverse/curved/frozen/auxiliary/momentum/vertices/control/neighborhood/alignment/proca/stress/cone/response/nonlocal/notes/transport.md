# Exact varied transport and full radial integrals

Use the unchanged selected mode's Bogoliubov coefficients relative to W8. The source is zero on an initial neighborhood, so their independent variations vanish initially. The coefficients themselves are not vacuum reset: the frozen preparation bound is

    |B_mix| <= C_mix/nu^6 < 1,
    C_mix=1813229+5347035781757616/m^4,
    nu²=m²+|p_com|²/Amax², Amax=25/16.

The analytic Wronskian gives |A_mix|<2. The exact mixing decomposition and its nonzero initial phase-interference control are checked independently.

Let C bound the old background residual, C_delta the newly recomputed varied residual, and B_freq the new bound on |delta W|/omega. Since I has length one, nu<=omega<=Amax nu and the transport exponential is below two, variation of constants gives

    |delta A_mix|,|delta B_mix| <= D9/nu^9 + D8/nu^8,
    D9=ceil[16(C_delta+2C B_freq)],
    D8=ceil[16C B_freq Amax].

The new common constants are D9=8296904784986635 and D8=260404337448359. They are recomputed from the new variations, not copied from the old mass-source response.

The reference readouts obey |J_ref|,|Z_ref|<=3 omega. Differentiating their moving canonical normalization and physical contacts before estimating gives a new common bound |delta J_ref|,|delta Z_ref|<=30 nu² per unit source norm. In particular the new physical contact envelopes are retained; no smaller old mass-only contact estimate is assumed.

The differentiated exact-minus-reference mixing identity is bounded by

    30 Amax nu(D9/nu^9+D8/nu^8) + 6 C_mix K_ref/nu^4.

The second term is essential despite zero independent initial response. The radial integrals of nu^-4, nu^-7 and nu^-8 are respectively Amax³/(8 pi m), Amax³/(15 pi² m^4) and Amax³/(64 pi m^5). The reference-tail radial integral is 1/(6 pi² m²). These identities are checked from the exact gamma-function integral, and all three physical polarizations are summed. No grid or ultraviolet cutoff replaces these integrals.

After dividing by L² and adding the new S6.84 physical finite local response, the complete vector energy and pressure bounds at L=10^400 are below 3.656e-795 and 4.697e-795 respectively. These decimal upper displays are rounded upward; the report uses exact rationals. The source norm is C10, including the C4 norm needed for the local part.

