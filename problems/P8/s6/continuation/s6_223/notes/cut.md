# Closed cut and endpoint bounds

Polynomial division of y^(2j)/(d-y^2), with d=1+4m^2/p, gives J0=atanh(1/sqrt(d))/sqrt(d) and Jj=d J_(j-1)-1/(2j-1). Therefore

A2(d)=-172/225+19d/30-d^2/10
       +sqrt(d)(30-20d+3d^2)atanh(1/sqrt(d))/30.

The formula is interpreted by analytic continuation from p>0; p0 is removable. On the upper bank p=-tau+i0, tau>4m^2, put z=1-4m^2/tau. The branch identity for atanh gives A2=D+i pi U, where

D(z)=-172/225+19z/30-z^2/10
     +sqrt(z)(30-20z+3z^2)atanh(sqrt(z))/30,
U(z)=sqrt(z)(30-20z+3z^2)/60.

U>0 on0<z<1, so neither bank can vanish even if D crosses zero. The reciprocal cut density is

rho(tau)=U(z)/[D(z)^2+pi^2 U(z)^2]>0.

At threshold, D0=-172/225 and U/sqrt(z)->1/2. D is analytic in z there because sqrt(z)atanh(sqrt(z)) is analytic. Thus rho=sqrt(z)h(z), with h analytic near0 and h0=50625/59168.

For a(Omega)=rho(Omega^2), Omega=2m+x, z=x/m+O(x^2/m^2). It follows, with differentiated remainders coming from this convergent local expansion, that

a(2m+x)=50625/(59168 sqrt(m)) sqrt(x)+O(x^(3/2)).
After subtracting its leading square-root term, the second derivative is O(x^(-1/2)) and is locally integrable.

At infinity use u=log(tau/(4m^2)), z=1-exp(-u), beta=sqrt(z). The exact stable identity

atanh(beta)=u/2+log(1+beta)

gives D as an explicit smooth function of u. Every u derivative of beta-1 and of its polynomial combinations is O(exp(-u)); the only growing factor is the explicit u. Consequently

D=(13/60)u+(13/30)log2-52/225+O(u exp(-u)),
D_u=13/60+O(u exp(-u)), D_uu=O(u exp(-u)),
U=13/60+O(exp(-u)), U_u,U_uu=O(exp(-u)).

These estimates follow by differentiating the explicit expression, not by differentiating an unsupported big-O assertion. With L=log(tau/m^2)=u+log4, the displayed constant becomes D=13L/60-52/225+O((m^2/tau)L). Hence rho L^2->60/13.

For a concrete tail lower bound, u>=16 implies z>=15/16, beta>=3/4 and30-20z+3z^2>=13. The polynomial part of D is greater than-1; log(1+beta)>=0. Therefore D>=13u/80-1>=u/10, while U<=1/2. Thus0<rho<=50/u^2 in this tail. The exact differentiated expression then gives rho_u=O(u^-3), rho_uu=O(u^-4). Changing from u to Omega yields

a=O((log Omega)^-2),
a'=O(1/[Omega(log Omega)^3]),
a''=O(1/[Omega^2(log Omega)^3]).

On compact interior intervals the denominator is nonzero because U>0. Together with the threshold expansion, these bounds give a(2m)=a(infinity)=0 and a' inL1(2m,infinity). They also give the second-derivative condition needed after subtracting the threshold model. No global numericalL1 constant or physical spectral cutoff is asserted.
