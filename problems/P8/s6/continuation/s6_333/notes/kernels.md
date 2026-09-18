# Remaining-energy kernels and endpoint control

Write q_e=(deltaa_e-deltaa)/(2e), g_e(y)=y^(2e),
h_e(y)=(g_e(y)-1)/(2e). Exact algebra gives
Z_e=q_e*g_e+deltaa*h_e. At e=0 the definitions are Delta,1,ln y.

For 0<y<=y+t<=x<1, 0<=g_e<=1 and |h_e(y)|<=|ln y|.
Integrating derivatives g_e'=2e*y^(2e-1) and h_e'=y^(2e-1)
gives the increment bounds 2e*ln(1+t/y) and ln(1+t/y).
They include the limiting e=0 cases.

For l=|ln y|, 1-exp(-2e*l)<=2e*l and
0<=exp(-2e*l)-1+2e*l<=(2e*l)^2/2.
Thus |g_e-1|<=2e*l and |h_e-ln y|<=e*l^2.
The latter Taylor inequality follows from zero value/first derivative
at zero and second derivative exp(-t)<=1, not from an uncontrolled
truncated series near y=0. The conditional moments in moments.md
make both endpoint logarithms integrable.
