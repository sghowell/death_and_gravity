# Complete nonforward finite part and specified soft prescription

Fix s strictly between4mu and n, and -1<z<1. Bounds below can be uniform
on compact subsets of this domain; no uniform forward limit is asserted.

Put d=w^2-1 and L=log(4/d). The exact master gives
J(w,w,z)=L/(1-z^2)+Fz+o(1),
Fz=[log((1-z)/2)/(1-z)+log((1+z)/2)/(1+z)]/2,
J(w,v,z)=L/[2(v^2-z^2)]+Hvz+o(1),
Hvz=[log((v-z)^2/(v^2-1))/(v-z)
    +log((v+z)^2/(v^2-1))/(v+z)]/(4v).
Also Q0(w)/w=L/2+o(1), Q2(w)/w=L/2-3/2+o(1).
For the equal pair write c=1-z>0, Delta=sqrt(c^2+2cd);
the logarithm is L+log[(c+d+Delta)/4]. Its smooth factor tends c/2.
For the mixed pair c tends v-z>0 and the logarithm is
L/2+log[(c+Delta)/(2sqrt(v^2-1))].
The complementary z->-z endpoints follow identically. Smooth inverse
prefactors introduce only O(d log d) errors on each compact domain.

The divergent coefficient is C0=P*A(z), and the complete finite part is

 Faux=P^2 Fz+2PB Hvz-3P a2 P2(z)+B^2 J(v,v,z)
      +2B[a0 Q0(v)/v+a2 Q2(v)P2(z)/v]+a0^2+a2^2 P2(z)/5.

To obtain a dimensional result, specify the minimal continuation of the
ACTUAL leading vacuum vertices to D=4+2epsilon, epsilon>0.
The conserved scalar-stress projector is trace-subtracted by1/(D-2).
Direct invariant contraction gives N_D=N4+2mu^2 epsilon/(1+epsilon)
in each gravitational channel. Therefore
A_D=A+epsilon A1+O(epsilon^2),
A1=-(2mu^2/kappa)(1/s+1/t+1/u),
P_D=P+epsilon P1+O(epsilon^2), P1=8mu^2/(kappa Q).
This is a regulator prescription for these vertices, not a unique
dimensional extension of the whole off-background DHOST action.

For the normalized S^(2+2epsilon) average, beta/gamma integration gives
<T_1(x)> = Gamma(3/2+epsilon)/(sqrt(pi)Gamma(1+epsilon))
           *B(1/2,epsilon) =1/(2epsilon)+1 EXACTLY.
Subtract P_D A_D(z)[T_1(a.n)+T_1(b.n)] from A_D(a.n)A_D(b.n).
For strictly nonforward z the four poles at n=+-a,+-b are disjoint.
Near each pole at angular distance theta, the other-axis tree differs
from its endpoint value by O(theta), uniformly for small epsilon.
The residual is bounded by a constant times1/theta plus a bounded
function. Its angular measure supplies theta^(1+2epsilon)dtheta.
This is integrable uniformly down to epsilon0, including after one
epsilon derivative with its logarithmic theta factor.

For noninteger dimensions this argument is made in two polar variables:
the marginal measures are normalized beta densities
(1-x^2)^epsilon dx and (1-y^2)^(epsilon-1/2)dy. On a neighborhood of
epsilon0 their weights and first derivatives admit the same integrable
majorants, after choosing a polar chart at each of the four poles.
Away from the poles the denominators are bounded. Dominated convergence
and its first-derivative bound therefore apply. For the auxiliary
w regulator, subtract P*A(z)[T_w(a.n)+T_w(b.n)] instead.
Its O(d)/(d+theta^2) endpoint mismatch vanishes in integral, and the
remaining residual has the same integrable limit. Hence that common
limit is exactly Faux, not an arbitrary finite-part choice.

Writing C1=P1*A+P*A1, the normalized dimensional average is
C0/epsilon+Faux+2C0+C1+O(epsilon).

Derive the two-particle measure from both on-shell delta functions:
dPhi2^D=(2pi)^(2-D) p^(D-3)/(4sqrt(s)) dOmega_(D-2).
Take each canonical four-point tree proportional to nu^(-2epsilon)
and rescale the final amplitude by nu^(2epsilon). Including both
halves, the cut prefactor is
beta/(32pi)*(s beta^2/(16pi nu^2))^epsilon
*Gamma(3/2)/Gamma(3/2+epsilon).
Using psi(3/2)=2-gamma_E-2log2 gives the whole finite bare cut

 rho_bare=beta/(32pi){C0/epsilon+Faux
          +C0[log(Q/(4pi nu^2))+gamma_E]+C1}+O(epsilon).

Now use precisely S278's DEFINED four-dimensional analytic soft kernel.
Its continuation satisfies
Im logW_E=beta P/(32pi)[1/epsilon+2log(E/nu)]+O(epsilon).
The same sign follows from the positive s+i0 discontinuity of S278's
pair parameter integral. Subtract this times A_D, not only A.
P*A1 cancels, but P1*A remains:

 rho_stripped=beta/(32pi){Faux
       +P*A[log(Q/(4pi E^2))+gamma_E]+P1*A}.

The dimensional scale cancels and the exact resolution derivative
matches S278. The finite gamma_E/log4pi and P1*A terms follow from the
explicit raw measure and minimal-D-versus-four-dimensional-soft choices.
A different finite soft convention must be converted explicitly;
neither an implicit MSbar rescaling nor identifying w with E is allowed.

The reference definition is Bellazzini et al.,2512.13780v2, equations
2.1--2.4 and3.8: https://arxiv.org/html/2512.13780v2 .
Its all-order full-factorization premise and finite-coupling limitations
are not proved here. This is a conditional first-loop cut in the stated
prescription, not a finite-G optical positivity theorem, a full complex
amplitude or a crossed forward subtraction.
