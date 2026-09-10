# Actual polynomial-vacuum light pole

Work in the canonical mass-one flat-space decoupled polynomial model,
not on the rolling background. Write M=mu^2 and g=G^2 for the heavy mass
squared and squared cubic coupling. All actual values are imported from
S6.110 and retain its fixed interactions and stable classical vacuum.

Define the Minkowski inverse propagator convention
D(s)=s-1+Pi(s). At one loop the only momentum-dependent contribution is
g/(16 pi^2) times the light-heavy Feynman-parameter bubble.
Choose Pi_R(s)=Pi(s)-Pi(1)-(s-1)Pi'(1) once.
The mass-one pole and unit residue refer to this perturbatively truncated
inverse, not to an exact all-orders propagator.

Prove its analytic parameter kernel, positive derivative, finite local
kinetic subtraction, conservative complex-disc remainder and zero exclusion.
Explicitly reconcile this choice with S6.110's zero-momentum potential
subtractions. Do not claim that pole mass and zero-momentum curvature can
both be set to one by the same mass counterterm.

The ordinary scientific calculation uses unmodified SymPy. Written
complex-analysis and continuum integral arguments are source-pinned, not
proof-assistant formalized. The exact checks are finite algebraic inputs
and rational bounds, not a machine proof of every continuum theorem.
