## Quantitative infrared cutoff of a and Delta

S330's relative conditioning loss is<=a*eta/x.
On the limiting cut, the positive-increment index difference has
mean<=1700a*eta/kappa. Changing the conditioning event adds
<=2*1400a*eta/kappa. Thus actual cutoff-conditioned index means
differ by<=4500a*eta/kappa, physically<4000eta/kappa^2.

For f(t)=t(1+ln(1/t)), subadditivity on positive energies gives
f(T_eta)<=sum_(w<eta)f(w) on the total-cut event.
The one-distinguished-emission identity bounds its conditional mean by
a*int_0^eta(1+ln(1/w))(1-w/x)^a dw
 <=a*eta*(2+ln(1/eta)).
The same-event Delta error is thus<=11000 times this quantity/kappa.
The event-change error is<=20000a*eta*(1+ln(1/x))/kappa.
Since eta<=x, the sum is<=42000a*eta*(1+ln(1/eta))/kappa,
physically<34000eta*(1+ln(1/eta))/kappa^2.

This does NOT estimate cutoff changes of the remaining-energy logarithm
in Z itself, does not promote state-dependent a to a branching rule,
and does not infer missing D-dimensional hard/evanescent matching.

For completeness, f(r)+f(s)-f(r+s)
=(r+s)ln(r+s)-rln(r)-sln(s)>=0 for r,s>0.
The finite subadditivity inequality passes to the infinite positive
tail by monotone sums on the total-cut event; the distinguished-emission
integral is finite. No Jensen replacement introduces ln(1/a).
Both compared state marks subtract the SAME Born value before the
nested-event mixture estimate. These are true finite-cutoff versus
limiting conditional means, not only a same-event comparison.
