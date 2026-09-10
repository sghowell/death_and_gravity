# Spectral mixed-mass bubble and complete outer on-shell subtraction

Use the S6.135 measure w(u), supported on u>=T=4mF^2, and g=G^2.
The mixed bubble's momentum-dependent inverse convention is
the frozen scalar Pi convention:

    Pi_spec(s)=(g/Q) integral_T^infinity w(u) B(M,u;s) du,
    Delta(s)=x u+(1-x)M-x(1-x)s.

The continued Euclidean inverse receives -Pi_spec. The prefactor
g/Q is checked against the frozen original H-Phi bubble.
Its mass-one light limit agrees after reflecting x -> 1-x.

The unsubtracted spectral expression has an overall local mass
divergence. Define the complete fixed on-shell difference before
performing that integral:

    Pi_spec,R(s)=Pi_spec(s)-Pi_spec(1)-(s-1)Pi_spec'(1).

For A=x(1-x), Delta_1=Delta(1), b=A/Delta_1, its finite parameter
kernel is

    R(s)=-log[1-(s-1)b]-(s-1)b.

The derivative and anchored value agree with the uncombined
logarithms, which identifies the analytic branch without a
global logarithm manipulation. The proper inner fermion
subtraction has already been paired in w and is not repeated.

For M>=1, u>=T>=16 and Re(s)<=4,

    Re Delta >= x(u-4)+(1-x)M>0.

The difference between the two sides is A(4-Re s)+4x^2.
This domain includes the entire radius-two disc about mass one.
The finite parameter representation follows first in real
Euclidean kinematics and then by the convergent analytic
continuation justified in the bounds note. No divergent
complex-momentum contour translation is assumed.

M is the same perturbative heavy mass reference. The calculation
does not declare H to be an exact stable asymptotic particle.
