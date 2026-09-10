# Full-Hessian normalization

The S6.110 exact Gaussian heavy integral has field-dependent light action
S_eff=Phi K_m Phi/2 + lambda4 sum(Phi_i^4)/24
      -G^2 J K_H^(-1) J/8, J_i=Phi_i^2.
Its second variation is K_m+W, where

W_ij = [lambda4 Phi_i^2/2
        -G^2 (K_H^(-1)J)_i/2] delta_ij
       -G^2 Phi_i (K_H^(-1))_ij Phi_j.

The new calculation differentiates the actual generic two-site effective
action before comparing all four entries with this formula. It independently
checks every contribution in Tr(K_m^(-1) W)/2. No mixed fluctuation determinant
is replaced by the field-independent heavy Gaussian determinant.

For translation-invariant continuum kernels, the first diagonal term is
the local quartic tadpole. The stationary-heavy term becomes
-G^2 G_m(0) integral(Phi^2)/(4M), also a momentum-independent light mass term.
A heavy one-point counterterm can redistribute this local term; it does not
alter the following nonlocal kernel or its on-shell subtraction. The mixed
term is -G^2 Tr(G_m Phi G_H Phi)/2. Thus the Euclidean inverse acquires
-G^2 times one light-heavy bubble. Continuing p_E^2=-s and changing from
the Euclidean to D(s)=s-1+Pi(s) convention gives

Pi(s)=g/(16 pi^2) [local UV constant - integral_0^1 log Delta(s,x) dx]
       + additional local affine terms,
Delta(s,x)=x M+1-x-x(1-x)s.

The factor follows from d^4q/(2pi)^4 = y dy/(16pi^2) after angular
integration, and the differentiated integral
integral_0^infinity y/(y+Delta)^3 dy=1/(2Delta).
The code verifies the literal four-component Feynman-parameter shift and
this convergent radial integral. It does not assign a finite value to an
unregulated divergent bubble. Differences and derivatives below are finite.

With Delta-i0 on the cut, -log(Delta-i0) has positive imaginary part.
This matches the stated inverse-propagator convention, not its opposite.
