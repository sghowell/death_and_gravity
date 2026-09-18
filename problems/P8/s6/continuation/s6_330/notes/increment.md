# Uniform positive-energy increment proof

Fix E,u and sigma=tau+rho positive, with R_sigma<=1/8 and t=massrho.
Follow Q(s)=Q_tau+sQ_rho. Positivity keeps this entire path in the
original physical recoil domain, with|Q_rho,sp|<=t.
Write K=(2E,0)-Q,m=sqrt(K^2),r=sqrt(m^2/4-1).
The original S300 outgoing momenta are exactly
p_pm^0=K0/2 +/- (r/m)(Ksp.u),
p_pm,sp=Ksp/2 +/- r[u+Ksp(Ksp.u)/(m(K0+m))].
This is a rewrite of the frozen recoil, not a new prescription.

Uniformly m>2,m<=4,r>3/5,r<2,|Ksp|<=1/8 and
|K0'|,|Ksp'|<=t. Since|m'|<3t and|r'|<5t, the boost
Bsp=Ksp(Ksp.u)/(m(K0+m)) obeys
|Bsp|<=1/512 and|Bsp'|<(1/32+5/512)t<t/16.
For the latter, the numerator derivative is<=t/4, the denominator
is>8 and its logarithmic derivative is<5t.
Thus|p_sp'|<t/2+5t*(1+1/512)+2t/16<6t.
Also|(r/m)'|<4t,r/m<1/2, so|p0'|<2t.
Incoming momenta do not change. Integration gives the same finite
difference constants.

For a unit spatial TT tensor A and q=(1,n), define D_i=k_i.q,
P_ij=A(k_i,k_j),L_i=P_ii/D_i and
S_ij=D_j L_i+D_i L_j-2P_ij.
The original compact domain has1/4<|D_i|<4,|L_i|<16,
massive spatial norms<2 and|S_MM|<102.
Outgoing differences obey|deltaD|<=8t,|deltaPii|<=24t,
|deltaLi|<=608t,|deltaS_MM|<=5168t and|delta z_MM|<=32t.
The last follows by differentiating each timelike dot product along
the physical path, not by leaving the mass shell.
With|f'|<=56,|f''|<161 on the frozen massive dot domain,
the six massive-pair contribution to deltaF is at most
6*(56*5168+102*161*32)t/144<34000t.

For an unweighted null ray v=(1,vhat), put d=v.q and ell=A(v,v)/d.
Off the soft-collinear point,|ell|<=2, and
s_i(v)=dL_i+D_i ell-2A(k_i,v) satisfies|s_i|<=44.
Only two outgoing massive legs vary;|deltas_i|<=1244t.
For alpha=|k_i.v|,|ln(2alpha)+1|<4 and|delta lnalpha|<=32t.
Existing mixed terms plus the newly added tail cost at most
[2*(4*1244+44*32)/8+4*4*44]t/36<65t.
The null-null Gram bound implies|[ln(2delta)+1]s|<=10.
The old-new plus new-new weighted mass is<=R_sigma*t, so the
null-null change costs< t. Relative-collinear limits are included.

For Sbar*H/(4pi^2), use|Sbar|<49,|H|<25,
|deltaSbar|<=2*608t+2t=1218t and|deltaH|<50t.
The massive part of H uses|d(Dln|D|)/dD|<5/2,
two varying massive legs and|deltaD|<=8t; the new null term
uses|d ln d|<2. Hence the real phase costs
<(1218*25+49*50)t/36<914t.
Together|deltaF|<35000t.

S328's aggregate phase is
G=[(c12+2)S12+(c34-2)S34]/(8pi).
Only the outgoing pair varies. Its gamma=(K^2-2)/2 has
|gamma'|<=33t/8<5t. On gamma>=29/16,
c'<2/5 and|c-2|<1/7. Therefore
|deltaG|<(102*2+5168/7)t/24<40t.
The complete norm is<35040t<40000t for t>0, and zero at t=0.

First work off all finite soft-collinear directions. S328 continuity
of the COMPLETE coefficient extends the bound to those directions;
no individual-current value is assigned. Approximate tau and rho
separately by positive atomic measures of the same masses. S328's
weak-to-uniform extension passes this addition-preserving estimate
to all positive Borel measures. It does not establish the same
bound for arbitrary signed perturbations or varying E,u.
