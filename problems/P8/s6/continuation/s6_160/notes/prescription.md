# Fixed map and matched physical source

At the named canonical MS interaction boundary, fix lambda, gamma
and c=4 lambda^2/gamma to the S6.111/S6.156-S6.159 numerical values.
There is no additional formal loop-parameter dependence in F.
Write E=Box+1, K=Box+2, X=(partial Psi)^2 and
Z=partial^mu Psi partial^nu Psi partial_mu partial_nu Psi. Then

    R = -c Psi^3/6 + lambda(2 Psi X + Psi^2 E Psi)
        + gamma(2 Z + Psi X + Psi^3/3)
        - gamma Psi K X - gamma Psi K(Psi E Psi)/2.

The frozen literal jet polynomial has 61 monomials, cubic field
degree and total derivative degree at most four. Its linear part
is the identity. Its inverse exists as a formal local field-degree
series; no global inverse or derivative-norm contraction follows.

Use D=4-2 epsilon, continue Lorentz contractions consistently and,
in light-mass-one units, multiply each of c, lambda and gamma by
mu^(2 epsilon), with mu=mF. This gives each term the D-dimensional
field dimension and preserves c_D=4 lambda_D^2/gamma_D. The
dimensionful constants in E and K keep their stated dimensions.
The same continuation is used everywhere before Laurent finite
parts; the four-dimensional jet enumeration is not a substitute
for the D-dimensional tensor prescription.

Define Z'[J,J_H] by substituting Phi=F_D(Psi) in the regulated
canonical parent functional, including every fixed parent
counterterm and det(delta F_D/delta Psi), with source
J F_D(Psi)+J_H H. Fermion and gauge integration variables are
unchanged. This defines the joint renormalized action/source
prescription by pullback. It is not independent minimal
subtraction of the ordinary Psi operator. Local composite
subgraphs can mix with fermion operators; this construction
does not claim their separate ordinary-Psi MS coefficients.

The reference values and F are held fixed under source
differentiation. If F were instead prescribed as a function of
running or bare couplings, its induced variations would have to
be included. That is a different dictionary, not silently used
here. There are no new freely adjusted b2 or vacuum parameters.
