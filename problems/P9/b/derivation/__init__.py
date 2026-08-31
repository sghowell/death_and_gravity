"""P9(b) S1 derivation chain (exact sympy; every step an identity test).

Frozen object: GLV14 (arXiv:1411.3712) eq. (86) five-alpha ADM quadratic action
(alpha_H: their eq. (82); alpha_M = d ln M^2 / dx).  Everything downstream
(background equations, D, c_s^2, I1, convention maps, dictionaries) is derived
here; published formulas BS14 (3.12)-(3.13), GLV14 (77), (79)-(80), (83), (85),
(88)-(89) are known-answer tests, not inputs.  House style: problems/P4/derivation.

Modules
-------
tools        symbols, eps-jet helpers, x-averaging, IBP canonicalizer
geometry     exact ADM geometry of the perturbed flat-FLRW scalar sector
actions      frozen (86) transcription; GPV-form (87); k-essence matter P(Y)
background   tadpole (eps^1) equations; family-B definitions; conservation
reduction    constraint elimination -> D, c_s^2; known-answer tests; I1; maps
kgb          covariant KGB in unitary gauge: independent route + DPSV pins
"""
