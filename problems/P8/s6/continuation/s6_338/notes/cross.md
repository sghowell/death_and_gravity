## Across the sources: retain EOM terms and both loop insertions

The complete cross-contraction is

2gc integral Phi(x)[Phi(y)Y(y) G_H G
       +Phi(y)^2 nabla Phi(y) G_H nabla_y G].

Both integrals carry the covariant volume elements. Use
Phi Y=Box(Phi^3)/6-Phi^2 Box Phi/2 and
Phi^2 nabla Phi=nabla(Phi^3)/3. Covariant integration by parts,
including every connection term, gives

gc/3 Phi(x)Phi(y)^3 [G Box G_H-G_H Box G]
  -gc Phi(x)Phi(y)^2 Box Phi(y) G_H G.

With (Box+n)G_H=-delta and (Box+1)G=-delta this is

gc(4-n)/3 Phi(x)Phi(y)^3 G_H G
  -gc/3 delta(G-G_H) Phi^4
  -gc Phi(x)Phi(y)^2 (Box+1)Phi(y) G_H G.

All 24 flat labels reproduce
32g^2/kappa [A0(1)-A0(n)+(n-4)B0(1;1,n)].
The last term is NOT set to zero before radiation: external emission
from its marked field produces -2 epsilon(p_i,p_i)F(1) for each
assignment, and varying Box produces the opposite contribution.
The complete 24-label external sum and metric sum are each nonzero
and cancel exactly. Varying other factors multiplies a vanishing
unshifted light EOM. No unproved field-redefinition shortcut is used.

The two massive propagators in K=G_H G each have a graviton insertion.
Write a=p.k, k^2=0, v=p^2, and
M(z;v)=1-z+n z-z(1-z)v.
Completing the square gives the light-line triangle denominator
M_L=M-2zya, y in [0,1-z], with TT weight z^2; and the heavy-line
denominator M_H=M-2(1-z)ya, y in [0,z], with weight (1-z)^2.
After radial integration, relative to Gamma(1+epsilon), their
integrands are -2z^2 M_L^(-1-epsilon) and
-2(1-z)^2 M_H^(-1-epsilon). Their exact antiderivatives are
-z M_L^(-epsilon)/(a epsilon) and
-(1-z) M_H^(-epsilon)/(a epsilon).
Summing endpoints and using Gamma(1+epsilon)=epsilon Gamma(epsilon)
gives the whole-D physical TT kernel insertion

-epsilon(p,p) [B_D((p+k)^2)-B_D(p^2)]/(p.k).

The pole of B_D is independent of p^2 and cancels in this difference.
The finite MS identity follows from B_MS(v)=-integral log M(z;v)dz.
This is an explicit two-triangle calculation, not an assumption that
a flat kernel has a unique arbitrary curved completion.

In the original recoil domain |p.k|<=4omega<=1/2, so the shifted
v lies in [0,2]. For n>=128 and 0<=z<=1,
M(z;v)>=1+(n-3)z+2z^2>=1; the triangle interpolants share this gap.
The kernel is analytic at null k. No massless 1/k^2 response is hidden.
A two-endpoint homogeneous symmetric response A pp+B(pk+kp)+Ckk+Deta
would have k contraction A(p.k)p+[B(p.k)+D]k; at p.k nonzero its
physical TT coefficient A vanishes. This check is specific to this
two-endpoint massive kernel, not a generic graviton-loop soft theorem.

For Phi^3 F(-Box)Phi the four external emissions are
6 sum_i J_i [F(1+2a_i)+3F(1)]. The kernel contact is
-6 sum_i J_i [F(1+2a_i)-F(1)]. Their total is
24F(1) sum_i J_i, precisely constant-quartic radiation.
The coincident tadpole contact also has only constant-quartic
physical TT response: its extra onepoint loop depends only on k,
and its linear TT response vanishes in exact D.
The fixed mixed UV counterfunctional Phi^2Y obeys the covariant
divergence/EOM identity Phi^2Y~Phi^4/3 for this complete amplitude;
all external and EOM metric emissions are retained as above.
Therefore the SAME existing OS4 linear constant subtraction removes
the entire across-source physical TT sector, not merely its flat value.
