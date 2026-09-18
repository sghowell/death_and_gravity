# Same-event remaining-energy logarithm

On A put y=x-R and y_eta=x-R_eta=y+T. The two original index bounds give

|deltaa_sigma*ln y-deltaa_eta*ln(y+T)|
 <=1700*T*|ln y|/kappa+1400*R_eta*ln(1+T/y)/kappa.

For 0<a<=1 Campbell's formula, applied to the unconditioned cloud,
and the exact conditional remaining-energy logarithm imply

E[T*|ln(x-R)| |A]
 =a*integral_0^eta (1-w/x)^a
   [-ln(x-w)+psi(a+1)+EulerGamma] dw.

The harmonic term is between 0 and 1. The integral of
-ln(1-w/x) from 0 to eta is <=eta: its ratio to x is
h+(1-h)*ln(1-h)<=h for h=eta/x. Dropping the nonnegative CDF ratio
therefore gives <=a*eta*(ln(1/x)+2)<=2*a*eta*L_eta.
The identity at eta=x agrees with the known R-weighted log moment.

For eta<=x/4, logarithmic subadditivity gives
ln(1+sum w/y)<=sum ln(1+w/y). Infinite sums are allowed: for each
y>0 their sum is <=T/y; nonnegative expectations use Tonelli and
monotone limits. Selecting one emission w, let q be all remaining
cloud energy. Its measure inside q<=x-w is dF(q)/F(x).
For q<=x/2 the log ratio ln[(x-q)/(x-w-q)] is <=4w/x.
For q>=x/2 its density is <=2a/x. The integral of the log ratio
in that region is <=I(w,x-w)<=w*[1+ln(x/w)], where
I(v,w)=(v+w)*ln(v+w)-v*ln v-w*ln w.
Integration against a*dw/w gives

E[ln(1+T/y)|A]<=a*eta/x*[4+2a*(2+ln(x/eta))].

After multiplication by 1400*R_eta<=1400*x, this is at most
11200*a*eta*L_eta/kappa.
For eta>x/4, use R_eta<=R and ln(y_eta/y)<=|ln y|, since
0<y<=y_eta<=x<1. The exact weighted moment is
a*x/(a+1)*[ln(1/x)+psi(a+2)+EulerGamma].
Its harmonic term is <=3/2, and x<4eta gives the smaller
8400*a*eta*L_eta/kappa. Thus both regimes fit 11200, and the
same-event total is <=14600*a*eta*L_eta/kappa.
