# Exact conditional logarithmic moments

For b=a+1, L_x=1+ln(1/x), l_x=ln(1/x), H_b=psi(b+1)+EulerGamma,

E[R*L_R^2|R<=x]
 =a*x/b*[L_x^2+2L_x/b+2/b^2] <=5a*x*L_x^2;

E[R*L_R*|ln(x-R)||R<=x]
 =a*x/b*[L_x*(l_x+H_b)+l_x/b+H_b/b-psi1(b+1)]
 <=3a*x*L_x^2;

E[R*|ln(x-R)|^2|R<=x]
 =a*x/b*[(l_x+H_b)^2+psi1(1)-psi1(b+1)]
 <=5a*x*L_x^2.

These are derivatives of the SAME conditional beta/moment integrals,
not new independent radiation laws. Use b>=1,H_b<=3/2 and
0<psi1(1)-psi1(b+1)<pi^2/6<2. For the last upper bound,
(l_x+3/2)^2+2<5L_x^2 when L_x>=1.
At a=0 all moments are defined as zero.
