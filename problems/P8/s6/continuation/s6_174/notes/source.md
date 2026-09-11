# Vacuum-regular lower dictionary and new clock alignment

Retain the unique analytic q from S6.109:
q_X+[1/(2X)-3R_X/(4R)]q=3R_X R_u/(4R).
Its coefficient-space integral is not a retarded solution in time.
In particular q(u,1) generally differs from zero. Resetting it to
the old clock boundary produces the previously proved vacuum
singularity; that old choice is not used.

The complete all64 stationary affine solve gives, in source-signature
-+++ notation, Tstar_mu=u_mu(alpha+c Box_source u+d uHu), with

    alpha=3(q+R_u/2)/2+3(R-1)R_u/(4R),
    c=-(R-1)/X,
    d=-(R-1)/X^2+3(R-1)R_X/(2RX).

These are obtained from the full source, including its coefficient
derivatives. The connection and P8 metric sign dictionary changes
Box_source to-Box_P8, while uHu is unchanged.

Take the NEW local coefficient

    B(u,X)=3H(u)(R-1)+3(q+R_u/2)/2,
    H(u)=4u/(1+u^2), W=T-B du, S=Tstar-B du.

This uses coefficient functions, not an external time source.
The q term cancels completely from the retained S:

    S_mu=u_mu(R-1)[-3H+3R_u/(4R)+Box_P8 u/X
                  +(-1/X^2+3R_X/(2RX))uHu].

R-1=O(X^2) removes every apparent X pole. The divided Hessian
coefficient generally has a NONZERO limit at a nonzero null
gradient; regular does not mean vanishing there. Around the
constant u=0 vacuum, S starts at four scalar fields. The new
vector-light mass-source vertex starts at total field degree five.

In the exact timelike metric chart from notes/chart.md,

    S_normal=(R-1)(K_hat-3H s), S_i=0.

No lapse velocity remains. Since R=1 on the entire clock,
K_hat=3H and s=1 there, S and every first variation vanish.
For s=1-epsilon*n and K_hat=3H+epsilon*k1,

    S_normal=epsilon^2[-2n(k1+3Hn)/h]+O(epsilon^3).

This is generally nonzero. The cubic interaction-W.S_second and
higher terms are retained; the vector is not a nonlinear spectator.

The full action therefore has the exact original classical clock,
with W=0, the old stationary connection, and original nonzero M1
momentum. Its metric/scalar equations agree with the complete
S6.173 target equations. The affine trace itself need not vanish:
on the clock Tstar=B du=3q(u,1)du/2. The shift is essential.
