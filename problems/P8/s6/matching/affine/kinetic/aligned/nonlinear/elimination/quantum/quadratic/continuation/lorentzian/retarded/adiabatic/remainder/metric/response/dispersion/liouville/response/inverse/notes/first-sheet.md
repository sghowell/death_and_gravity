# First-sheet inverse and positive cut density

All matrices are in the physical source order (N,Z). Fix h,m>0.
Set alpha=4/(9h), beta=28/(81h), c=16/(81h),

bT=(beta(1-z),2(1-z))^T,
bL=(beta+alpha*z,2(1+z))^T,
M(z)=2*bT*bT^T+bL*bL^T,
Q=[[1,0],[-c,1]],
F=[[1256/(6561h^2),-76/(243h)],[-76/(243h),-4]].

The exact integral in FORMULATION.md is holomorphic on
C minus (-infinity,-4m^2]. At p=0 it equals F without a pole.
Its integrable endpoint at p=-4m^2 is considered below.
The eigenvalue arguments concern complex vectors with Hermitian
inner product, not only real test vectors.

## No omitted zeros

For 0<z<1, det M=512*z^2*(1-z)^2/(6561h^2)>0 and M is a Gram
matrix. Therefore it is positive definite. Direct scalar
division gives, for nonreal p,

Im B(p)=-m^2 Im(p) integral_0^1
 y^2 M(y^2)/abs(4m^2+(1-y^2)p)^2 dy.

The matrix multiplying -Im(p) is strictly positive definite.
If B(p)v=0 for a nonzero complex v, its Hermitian imaginary
quadratic form would vanish, a contradiction. Thus all nonreal
first-sheet p are zero-free.

For p>=0, the S6.69 Q chart has A>=a*=15616/(98415h^2)>0 and
C<=-4, hence det B<0. This uses the exact massive integral, not
the leading logarithm.

For -4m^2<=p<=0 put x=-p/(4m^2). The added matrix is

B(p)-F=(1/4) integral_0^1
 [x*y^2/(1-x+x*y^2)] M(y^2) dy.

The scalar weight lies in [0,1]; the endpoint x=1 is interpreted
by continuity at the measure-zero point y=0. Consequently
F<=B(p)<=T in Loewner order, where

T=F+(1/4) integral_0^1 M(y^2)dy
 =[[28408/(98415h^2),268/(1215h)],
   [268/(1215h),-16/15]].

In the Q chart A>=3128/(19683h^2)>0 and C<=-16/15<0.
So this entire real subthreshold interval is zero-free as well.
In particular det T=-526352/(1476225h^2) is nonzero. Continuity
at the complex threshold follows also from the elementary
radial-moment formula, whose only new term is O(sqrt(p+4m^2)).

These arguments exhaust the first sheet. Since s->p=s^2 maps
Re(s)>0 to the slit plane excluding the negative real axis,
R(s^2)=B(s^2)^-1 has no right-half-plane pole.

## Both cut banks, including the sign

At p=-Omega^2+i0, Omega>2m, put z=1-4m^2/Omega^2 in (0,1).
Write M(z)=M0+M1*z+M2*z^2. The elementary S6.69 moments give

B(-Omega^2+i0)=D(z)-i*pi*W(z),
W(z)=sqrt(z)*M(z)/8,
D(z)=F+[M0+M1*(z+1/3)+M2*(z^2+z/3+1/5)]/4
       -sqrt(z)*M(z)*atanh(sqrt(z))/4.

This sign can also be checked before integration:
Im[-p/(4*(4m^2+(1-y^2)p))] is negative on the upper bank.
The open-cut weight is positive definite, with
det W=8*z^3*(1-z)^2/(6561h^2)>0. Hence both cut banks are
invertible; the lower bank is the conjugate upper bank.

Let R=B^-1 and define rho(tau)=Im R(-tau+i0)/pi. The matrix
inverse identity gives rho=R^dagger W R, strictly positive
definite on tau>4m^2. For a real symmetric two-by-two pair D,W,

X=det D-pi^2 det W,
Y=D11*W22+D22*W11-2*D12*W12,
rho=[Y*adj D-X*adj W]/(X^2+pi^2*Y^2)
    =[adj D W adj D+pi^2 det W adj W]/(X^2+pi^2*Y^2).

The denominator is abs(det B)^2 and is never zero on this open
cut. The second form explicitly displays positivity. The
certificate checks these generic identities independently of
the particular mass vertices.

## Infinity and the Cauchy representation

The exact elementary moments, uniformly in the exterior of the
first-sheet slit plane, imply

B(p)=E-P*log(p/m^2)+O(m^2/p*log(p/m^2)),
P=2*(c,1)^T*(c,1), E=F+C_large.

This is now used at infinity only. In the Q chart,
Q^T E Q=[[a*,b*],[b*,14/15]], b*=604/(1215h).
The logarithm acts only on the second chart coordinate.
Schur inversion gives

Rinf=Q*diag(1/a*,0)*Q^T
 =[[98415h^2/15616,-1215h/976],
   [-1215h/976,15/61]],
R(p)-Rinf=O(1/log(abs(p))).

The uniform asymptotic follows directly from the finite
polynomial radial formulas: their only non-polynomial term is
the principal logarithm, and the residual is O(p^-1 log p),
also on either exterior cut bank. Since a*>0, the Schur
denominator has magnitude comparable to log(abs(p)) there.

R-Rinf is holomorphic with no missing pole residues. A keyhole
Cauchy contour has vanishing outer contribution O(1/log R);
the threshold circle vanishes since R is bounded there.
The jump is 2*pi*i*rho. Changing the cut coordinate x=-tau
fixes the minus sign:

R(p)=Rinf-integral_(4m^2)^infinity rho(tau)/(p+tau) dtau.

The integral is absolutely convergent off the cut:
rho=O(sqrt(tau-4m^2)) at threshold and O(log(tau)^-2) at infinity.
It extends to both banks in the usual boundary-value sense.
There is no separately inserted subtraction constant or
undetermined residue. At p=0 it yields the exact sum rule

integral_(4m^2)^infinity rho(tau)/tau dtau = Rinf-F^-1 = Mmoment.

Mmoment has positive first minor and determinant
177147*h^2/99536384. Its complete rational entries are recorded
in the certificate.
