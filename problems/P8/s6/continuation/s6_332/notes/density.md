# Common-space density split and finite conversion

Since A is a subset of A_eta and p/p_eta=1-r, on A
M_eta/p_eta-M/p=((M_eta-M)-r*M_eta)/p.
After absolute values and expectation, adding the excluded event gives

E|1_A_eta*M_eta/p_eta-1_A*M/p|
 <=E[|M_eta-M| |A]+r*E[|M_eta| |A]
   +E[|M_eta|*1_(A_eta minus A)|A_eta].

Both event-normalization errors are necessary. For a generic constant
mark and a proper nested event, each contributes r and the L1 distance
is 2r. A difference of means would entirely miss this example. The
constant-mark example is an algebraic negative control, not a physical
Born-subtracted soft current.

For Delta, the same-event term is at most
11000*E[T*(1+ln(1/T))|A]/kappa.
Writing h(t)=t*(1+ln(1/t)), entropy subadditivity gives
h(sum w)<=sum h(w). Campbell's distinguished-emission formula therefore
bounds the conditional expectation by a*eta*(2+ln(1/eta)); its CDF
ratio is <=1. This does not assume independence after conditioning.

On either applicable cut, the retained state's Born difference is
<=10000*x*(1+ln(1/x))/kappa. Both event terms are bounded using
r<=a*eta/x. The full Delta density error is thus <=42000*a*eta*L_eta/kappa.

The log same-event budget from tail.md is 14600*a*eta*L_eta/kappa.
The log event-change budget from endpoint.md is
6300*a*eta*L_eta^2/kappa. Since L_eta>=1,
42000+14600+6300=62900<63000 gives the result. Multiplication by 4/5
gives 50400<51000 for the original physical bound. No x^a expansion
or division by an uncontrolled small absolute cut probability occurs.
