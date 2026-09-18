# Relative rare-event normalization and marked convergence rates

The omitted tail T_eta and the retained state are independent BEFORE
conditioning. This gives
P(R<=x)/P_eta(x)=E[P_eta(x-T_eta)/P_eta(x)].
Integrating dlnP_eta/dr<=a/r between x-T_eta and x yields
P_eta(x-T_eta)/P_eta(x)>=(1-T_eta/x)_+^a for0<a<=1.
For T_eta=x the right side is zero, so the zero-energy atom causes
no exception. For T_eta>x both the extended CDF and lower bound are
zero. Since(1-z)_+^a>=1-z and E[T_eta]=a*eta,
0<=1-P(R<=x)/P_eta(x)<=a*eta/x.
At a=0 both probabilities are one. This is a RELATIVE estimate
that remains valid for arbitrarily small unexpanded F(a)x^a.

Let A={R<=x}, A_eta={R_eta<=x}, with A subset A_eta.
Write D=C-C_Born and B=23667. On their appropriate events,
||D_sigma||sup,||D_sigma_eta||sup<=B*x.
The conditional law on A_eta is a mixture of its restrictions to
A and A_eta minus A, with excluded weight
r=1-P(A)/P(A_eta)<=a*eta/x.
Subtracting the same Born mark first bounds the difference of
finite-state conditional means on these events by2B*x*r.
On A, positive-increment control and the exact missing-energy moment
bound the further full-state difference by40000*a*eta.
Thus the uniform angular TT norm of the complete conditioned mean
difference is<=(40000+2B)*a*eta=87334*a*eta.

For squared uniform norms about Born, use
| ||D_sigma||sup^2-||D_sigma_eta||sup^2 |
 <=2B*x*||D_sigma-D_sigma_eta||sup on A.
The same-event expectation is<=2B*x*40000*a*eta.
The mixture variable is nonnegative and bounded by B^2*x^2,
so its event-change error is<=B^2*x^2*r, with no extra factor2.
Total error is<=(2B*40000+B^2)*a*x*eta.
The argument also applies at each angular/polarization point.

With physical normalization kappa^(-3/2) and original a<4/(5kappa),
the two errors are<70000*eta/kappa^(5/2) and
<2000000000*x*eta/kappa^4. They hold uniformly for
0<eta<=x<=1/8, preserve the complete complex phase and do not
expand x^a. These are marked D4 leading-reference estimates,
not D-dimensional evanescent matching or interacting dynamics.
