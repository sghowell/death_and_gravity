# Exact missing-energy moments and the renewal identity

All Poisson statements refer to the unconditioned fixed-Born
intensity Lambda(dn)*dw/w, with total angular mass a.
The conditioned process is not assumed to remain Poisson.

Distinguish one energy-weighted emission in the positive finite-cutoff
Bose series: n/n!=1/(n-1)!. Passing to the monotone low-energy limit
gives E[T_eta*1_(R<=x)]=a*int_0^eta P(R<=x-w)dw.
Using the exact CDF F(a)r^a for0<r<=1 and eta<=x,
E[T_eta|R<=x]=a*x/(a+1)*[1-(1-eta/x)^(a+1)]<=a*eta.
Both numerator and denominator contain the SAME F(a)x^a; no
absolute error is divided by a tiny probability.

Distinguishing two DIFFERENT emissions uses n(n-1)/n!=1/(n-2)!.
The square also has its diagonal, so
E[T_eta^2|R<=x]
 =a*int_0^eta w(1-w/x)^a dw
 +a^2*int_0^eta int_0^eta (1-(w+v)/x)_+^a dw dv.
With h=eta/x, this equals
a*x^2*[(1-(1-h)^(a+1))/(a+1)
       -(1-(1-h)^(a+2))/(a+2)]
+a^2*x^2*[1-2(1-h)^(a+2)+(1-2h)_+^(a+2)]/((a+1)(a+2)).
The positive-part term is essential when2eta>x.
The upper bound is a*eta^2/2+a^2*eta^2.
At eta=x these become a*x/(a+1) and a*x^2/(a+2),
not the square of the mean. At a=0 both are exactly zero.
The positive increment theorem consequently gives
E[||C_sigma-C_sigma_eta||sup|R<=x]<=40000*a*eta.

For the finite-cutoff total energy, P_eta(r)=P(R_eta<=r)
has zero-count atom eta^a at zero. Distinguishing an
energy-weighted emission in the finite Bose density gives
r P_eta'(r)=a P_eta(r-eta) almost everywhere on0<r<1,
where P_eta(s)=0 for s<0. There is no upper-end term because
the reference upper energy is1 and r<1. In particular P_eta
is constant on0<r<eta and dlnP_eta/dr<=a/r elsewhere.
The atom is included in this identity and in every conditioning
normalization. At a=0 handle the deterministic zero state directly.
