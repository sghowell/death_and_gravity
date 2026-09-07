# Why the higher-operator complement must be fixed

Let v_mu=partial_mu phi, X=v_mu v^mu, T=Box(phi),
V=v^mu phi_munu v^nu and Rv=v^mu partial_mu R_B. The product rule and
partial_mu X=2 phi_munu v^nu give the exact identity

    div(W R_B X v)
      = W_phi R_B X^2 + W X Rv + W R_B X T + 2W R_B V.

This is a boundary for compactly supported variations. If its first
term alone is assigned to F2, its apparent shifts are
Delta F2=W_phi X^2, Delta A1=0 and Delta Xi=-4W_phi X. The remaining
three terms are part of the higher-operator complement, not zero.
Thus a projected Xi is meaningful only together with the specified
representative and complement. This does not challenge the declared
retained-order basis or the independent tensor-symbol test of S6.30.

The literal physical tensor metric diag(1,-exp(q(t)),-exp(-q(t)),-1)
has unit volume and R_B=-q'^2/2. For phi=t and W=t the boundary density
is R_B+t R_B'=d(t R_B)/dt. The first term has Euler expression q'',
while the second has -q'', despite its coefficient t vanishing at the
center. The full Euler expression is identically zero. Evaluating
coefficient values at the center before varying fails this exact test.

This issue is present at the natural scales of the actual S6.20 parent.
In theta=T units with tau=1, set a=MF2*r/beta1, ell=(log r)_theta and
W=MF2*r^2*a*a_theta*ell_theta. At theta=0,

    r=2, a=c(c-2)/16, a_theta=ell=0,
    a_thetatheta=(13c^2-22c+8)/8,
    ell_theta=-4(c+2)/c,
    W=0,
    W_theta=-MF2(c-2)(c+2)(13c^2-22c+8)/8.

Consequently W_theta/(c-2) tends to -8MF2 as c tends to 2 from above.
A center-only order-(c-2)^2 potential projection need not retain that
scaling after redistributing these operators by a boundary. This is
not a calculated physical remainder, an action at c=2, or a full
degree-six Xi. The separately tracked full Euler operator is the
relevant check on a given tensor probe.
