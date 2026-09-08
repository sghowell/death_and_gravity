# Complete quadratic constraints and the physical matter cone

Use the same normalized variables as S6.41, signature -+++, q>0.
Let n, b and t be the lapse, scalar shift and original temporal trace;
sigma is the longitudinal trace potential. Before any field shift the
literal quadratic density is

    L/a³=L_CD_base+(t+d*n)²/2-q*sigma²/2
                    +zeta*q*(sigma_dot-t-d*n)²/2.

The mass source and the electric source have the same sign. The
independent tests solve all three original auxiliary equations directly.
Equivalently, set t=t_new-d*n, a unit-Jacobian point transformation.
Then the full density separates as

    L/a³=L_CD_base+t_new²/2-q*sigma²/2
                         +zeta*q*(sigma_dot-t_new)²/2.

All coefficient and source variations are retained in deriving this
quadratic action: Sstar starts at order two, while W starts at order
one. Variations of its quadratic mass tensor or the metric first
contribute at cubic order. The Maxwell term is minimally coupled to
the original physical metric in W variables; its background vanishes.
Thus no metric or free-matter quadratic term has been discarded.

## Constraints, positivity and crossing

The joint auxiliary determinant is
-4q²*Theta²*(1+zeta*q). Away from Theta=0, the temporal solution is
t_new=zeta*q*sigma_dot/(1+zeta*q), and the shift equation gives the
unchanged n=(v_dot+ell*s/2)/Theta. The vector kinetic coefficient is
C=zeta*q/[2(1+zeta*q)]>0. The complete kinetic matrix is diag(K_CD,C),
with the old positive square decomposition

    (s_dot+w*v_dot/Theta)²/2+J*v_dot²/Theta²+C*sigma_dot².

Here s denotes the free-matter perturbation, not sqrt(-x) from the
nonlinear source note. For a regular first-order description at the
crossing, keep the old canonical momenta and let pi be the vector
momentum. Exactly

    H=H_CD+(1+1/(zeta*q))*pi²/2+q*sigma²/2,
    H_CD=C0+R0²/J, n=-R0/J.

No inverse Theta occurs. The complete old physical scalar system and
its regular continuation therefore remain unchanged; the new canonical
vector block is regular for all finite time. At u=0 the old q>6
velocity chart may still be used, but it is not imposed on the regular
first-order system. The zero-momentum vector is handled separately.

## Original cone condition, not a substitute

The old two-scalar matrices obey K0=G0 exactly after all constraints.
The new longitudinal block has C->1/2 at high momentum and principal
gradient coefficient 1/2. Consequently

    K_infinity=diag(K0,1/2)=G_infinity,
    det(G_infinity-c²*K_infinity)/det(K_infinity)=(1-c²)^3.

The two transverse Maxwell polarizations also have principal speed
one in the original matter metric. The tensor block is unchanged.
This directly passes the original K>0, G>0 and K-G>=0 cone test at
quadratic rolling order, unlike the unshifted S6.41 curl. No assertion
about a massive finite-q phase velocity is used. Nor are the nonzero
curl principal limit and the zero-curl auxiliary rank limit exchanged.

This does not establish a nonlinear secondary-constraint rank, an
interacting cutoff, a vacuum S-matrix, or any adopted V/G/B UV gate.
