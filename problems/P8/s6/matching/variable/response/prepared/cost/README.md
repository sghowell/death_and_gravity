# S6.29: actual physical preparation cost

This child bounds what it costs to prepare the exact regular-light tensor
sector of [S6.23](../README.md) using the unchanged full parent's actual
g-metric source. It supplies neither a vacuum nor a local EFT.

On J=(-1/50,-1/100), K in [1,4], with all four initial physical data
zero, the full loading map C_K has the continuously certified Gramian

    9/(4*10^12) I <= C_K C_K* <=16/25 I.

The target is in the full Frobenius-normalized endpoint frame. Raw
canonical outer data and physical g data are not substituted for it.
An explicit degree11 physical-field Hermite construction gives a
zero-extended H0^2 source. For a nonzero full target z its sufficient
bounds are ||sigma||2<2*10^8||z|| and
||sigma_uu||2<4*10^14||z||. Exact light symplectic moments give the
necessary costs ||sigma||2>19||v|| for nonzero regular targets Lv,
and >900 for the unit even target.

With Fourier convention sigmahat(omega)=integral exp(-i omega u)sigma(u)du,
every such compact source has at least two thirds of its squared norm
outside any band |omega|<=Omega<=100. For unit v and source budget10^6,
the exact spectral optimizer in L2 has strict budget slack. The H0^2
infimum equals that L2 minimum and lies strictly above722/3 and at most
4*10^12/9 (above540000 for the even target). Attainment in H0^2 is not
assumed. The optimizer is characterized analytically, not numerically
computed. Nine actual loading moments needed for a sharper band-Gramian
calculation remain uncomputed.

The same fixed source at positive delta has normalized output error

    <=42e+delta(8600||v||+12600000 S)

relative to Lv, where e is its limiting loading error and S its actual
L2 source budget. The exact analytic prepared target instead has the
coefficient8400. The unit-target example S=10^6, e=0, delta=10^-17
has error below1/1000. This is not a bound on raw physical derivatives
and is not small throughout the whole delta<=10^-9 box.

See [FORMULATION.md](FORMULATION.md), the
[validated Gramian proof](notes/gramian.md),
[actual-source construction](notes/construction.md),
[optimization and Fourier proof](notes/optimization.md), and
[source/scope record](notes/sources.md). The certificate is
[physical-preparation-cost.json](certificates/physical-preparation-cost.json).

A fresh ordinary replay, from the repository root:

    PYTHONHASHSEED=0 .venv/bin/python -c '
    import runpy,sys
    from pathlib import Path
    from sympy.core.random import seed
    seed(0)
    for path in Path("problems/P8").rglob("src"):
        sys.path.insert(0,str(path))
    sys.argv=["p8_preparation_cost.verify","--check"]
    runpy.run_module("p8_preparation_cost.verify",run_name="__main__")
    '

For ordinary tests use the same two RNG seeds and run pytest on this
child's tests with -p no:faulthandler. The broad P8 exact-GCD adapter is
not used for ordinary certificate verification.

The source is sigma=tau^2 Pi/M^2, with Pi an infinitesimal conserved
external TT anisotropic stress. Its squared norm is not positive matter
energy. Source concentration is not a rolling spectral gap, physical
cutoff, nonlinear backreaction estimate, arbitrary-source light-only
reduction, original C/D operator dictionary or a UV verdict.
Original S6/P8 remain open.
