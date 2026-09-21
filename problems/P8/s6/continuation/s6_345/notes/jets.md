# Literal metric response and the third scalar jet

For a Gram word the algebraic inverse-metric variation replaces each
G_ij by-2*H_ij, where H_ij=epsilon(p_i,p_j).
Connection variations must also be retained.

A Hessian on field i with slots contracted into p_j,p_k contributes

    -[a_j H_ik+a_k H_ij-a_i H_jk],

times its remaining edge, with a_i=k.p_i.
A self-traced Hessian has contracted connection
2*k_mu*epsilon^(mu nu)-k^nu*tr(epsilon), which vanishes only
for the stated real TT graviton.

For the symmetrized third scalar jet, directly expand
partialGamma*partialPhi and the three Gamma*partialpartialPhi terms,
then symmetrize all three indices. All64 spacetime components are
checked against

    delta J3(u)=-(3*p.u+k.u)
       [2*k.u*epsilon(p,u)-(k.p)*epsilon(u,u)],

with the Fourier factors stripped consistently with the flat vertex.
Thus the triple-edge contact is

    -6Gij^2 Hij
    -(3Gij+a_j)(2a_j Hij-a_i Hjj)
    -(3Gij+a_i)(2a_i Hij-a_j Hii).

The k-dependent piece is required; pure metric differentiation or
ignoring the symmetrized connection would change the curvature answer.

The generic second-jet implementation is also compared with the whole
existing S337 Galileon contact, including its metric and connection
terms, before scalar shell specialization. This checks the original
relative vertex convention independently of the new box calculation.
All24 scalar label assignments are taken for every basis word.

