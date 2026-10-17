"""P2: certified resonance-free strip for subextremal RNdS (stage S1 onwards).

Modules (FORMULATION.md section 4):
  arbutil      precision context and small Arb helpers
  params       exact parameter points P_A, P_B -> certified Arb balls (params.json)
  p20_identity exact sympy identity chain for P2-0 (closed upper half-plane) and D4
  qnm_float    float (mpmath) QNM solver: Frobenius series + Taylor continuation
  spectral     independent float solver: hyperboloidal Chebyshev collocation (KAT-7)
  arbw         validated Wronskian W_l(omega) in Arb (Frobenius/Taylor with tail bounds)
  winding      validated small-box winding numbers and Krawczyk enclosures (KAT-6)
  kats         drivers for the known-answer tests KAT-1 ... KAT-7
  lg_bounds    P2-2: Liouville-Green thresholds W_l with Arb error-control integrals
"""
