# Complete analytic source on canonical jet polydiscs

Use the exact S6.177 canonical family, with kappa>=kappa0=10^800,
N=1024 and m=1000. The new S6.182 fixed F correction does not
change R, its affine source or the vector operator.
Let Y=eta^(mu nu)Phi_mu Phi_nu and Z=Phi^mu Phi^nu Phi_mu nu.
The complete source is S_mu=U_mu/kappa, with

    U_mu=V_mu+B_mu/sqrt(kappa)+C_mu/(kappa R_k),
    V_mu=Phi_mu a(Y BoxPhi-Z),
    B_mu=-3Phi_mu Y^2 a Href,
    C_mu=3Phi_mu Y^2 a[
           Y^2 a_Phi/4+(2a+Y a_Y)Z/2],
    a(Phi,Y)=[(R_base(u,X)-1)/X^2]/kappa0,
    u=Phi/sqrt(kappa0), X=Y/kappa0,
    R_k=1+Y^2 a/kappa, Href=4u/(1+u^2).

The all-Y regular formula is checked against the literal frozen
source, not its first vacuum Taylor term. Neither X nor Y is
assumed positive, and the null-gradient limit is retained.

For an outer complex domain |u|<=1/4, |X|<=rho=1/(16N),
write T=X^N/[X^N+(1-X)^N] and w=N X^2 exp(-N X^2).
The denominator modulus is at least
(1-rho)^N-rho^N>=15/16-rho^N>7/8.
Thus |T|<=2|X|^N<=|X|^2.
Also |exp(-N X^2)|<=exp(N rho^2)<2, since
N rho^2<=1/1024 and exp(x)<=1/(1-x) for 0<=x<1.
For B=T+(1-T)w, |B|<=5N|X|^2.
Consequently the removable analytic coefficient
A=(R_base-1)/X^2=B(X)(X-1)/[X^2(1+u^2)^3]
obeys |A|<=20N. This includes X=0 by analyticity.

On the inner domain |u|<=1/8, |X|<=rho/2, Cauchy gives
|A_u|<=160N and |A_X|<=640N^2.
Consider all15 independent scalar jet coordinates (one value,
four first derivatives, ten symmetric second derivatives),
each with complex modulus<=2. Then |Y|<=16, |BoxPhi|<=8,
|Z|<=128, |Phi_mu|<=2. The actual hierarchy kappa0>=512N
places their u,X images in that inner domain. Thus

    |a|<=20N/kappa0,
    |a_Phi|<=160N/kappa0^(3/2),
    |a_Y|<=640N^2/kappa0^2, |Href|<=10/sqrt(kappa0).

The complex R_k satisfies |R_k-1|<=5120N/kappa0^2<1/2.
The three complete U terms are bounded per component by

    10240N/kappa0,
    307200N/kappa0^2,
    943718400N^2/kappa0^3.

For the last term, 2|a|+16|a_Y|<=80N/kappa0, and the
bracket in C is at most15360N/kappa0. These inequalities
use kappa>=kappa0>=512N; they include all denominator terms.
After dividing their sum by N/kappa0, the maximum at the
minimal hierarchy is10240+4200/N<=11290<20000.

Therefore every S component has modulus below
M=20000N/(kappa*kappa0) on the full complex jet polydisc.
The unit polydisc about every REAL jet with components<=1
fits inside it. Cauchy bounds first jet derivatives by M,
and every second jet derivative by2M. No higher field
terms or generated Fourier harmonics have been omitted.
