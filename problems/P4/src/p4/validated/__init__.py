"""Validated (Arb ball-arithmetic) local expansions for the P4 CSS problem (Stage S2).

Modules
    arbseries   truncated power series over ``flint.arb`` (wrapping ``arb_poly``)
    polysys     polynomial systems  P(t,u) * (theta u) = Q(t,u)  with exact fmpq_mpoly data
    recursion   certified order-by-order solution of the coefficient recursion
    tailbound   rigorous geometric tail bound (Banach fixed point in l^1_nu on the tail)
    sonic       sonic-point expansion (Theorem A, item A1)
    centre      regular-centre expansion in powers of e^x (Theorem A, item A2)
"""
