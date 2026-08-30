"""Modelling-scope certificate for P4 Theorem A: the KHA CSS system *is* the
continuously-self-similar reduction of spherically symmetric Einstein-Euler, p = rho/3.

Modules
    einstein_euler     from-scratch sympy derivation (metric -> Christoffel -> Ricci ->
                       Einstein; perfect fluid; div T = 0; CSS ansatz; KHA variables)
                       and exact-identity comparison with ``p4.css.coeffs`` and with the
                       certified polynomial system ``p4.validated.systems.sonic_system``
    bianchi            the angular Einstein equation is implied by the others (exact)
    qjet               exact second-order jets over Q (forward-mode AD with Fractions)
    independent_check  the same reduction re-verified pointwise in exact rational
                       arithmetic, sharing no code with ``einstein_euler``
See notes/modelling-scope.md.  Run: uv run --with sympy pytest problems/P4/tests/test_derivation.py
"""
