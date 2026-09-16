# Physical normal-sheet contour and complete radial moments

For25/4<=s<=16 put v=x+i f(x), f=x(1-x)(1-2x),0<=x<=1.
The denominator is

A=1-s[x(1-x)+f^2]-i s x(1-x)(1-2x)^2.

Replacing f by eta*f,0<eta<=1, keeps the imaginary part negative away
from the endpoints and center. At the center A=1-s/4 is nonzero.
Consequently this deformation passes above the first Feynman root and
below the second, without crossing a zero, and implements A-i0.
At the negative-real contact use the continuous lower-bank logarithm.

Let r=(2x-1)^2 and q=1-r. Then
Re A=1-sq/4-srq^2/16, Im A=-srq/4.
For r<=1/8, Re A<=-47/128. For r>=7/8, Re A>=31/64.
On the middle interval r(1-r)>=7/64, hence|Im A|>=175/1024.
Thus |A|>1/8. Also Re[v(1-v)]<=5/16 and
|Im[v(1-v)]|<=1/4, whose squared bounds sum to41/256<1/4.
So |v(1-v)|<1/2 and |A|<9<10. Since f' lies in[-1/2,1],
the path length is below2. The logarithm has magnitude below8 using
ln8,ln10<3 and pi<4. Therefore

|J0|<16, |J1|<128, |Lll|<16, |Bbar|<8,

where Bbar=-integral v(1-v)ln A. For spacelike b in[-12,0],
1<=A_b<=4 on the real path and the same bounds are immediate.
These are uniform functional bounds, not a finite-grid claim.

## Heavy radial integral without an unresolved physical pole

For every A on that contour and n>=128 factor

n(1-u)+Au^2=n(1-alpha u)(1-beta u),
alpha=(1+sqrt(1-4A/n))/2, beta=(A/n)/alpha.

The real part of the square root exceeds4/5. Hence|beta|<1/10,
9/10<|alpha|<11/10, and|alpha-beta|>4/5. beta has the same lower
bank as A, and |ln beta|<ln n+7. For k=0,1,2 define

Q_k(c)=[-ln(1-c)-sum_(l=1)^k c^l/l]/c^(k+1),
Q_k(0)=1/(k+1).

The exact partial fractions give
I_j(A)=integral_0^1 u^j/[n(1-u)+Au^2]du
 =[Q_(j-1)(alpha)-Q_(j-1)(beta)]/[n(alpha-beta)], j=1,2,3.

Use|Q_k(beta)|<=10/9, (10/9)^3<7/5 and
11/10+(11/10)^2/2<2. Then

|I_j|<(5/4)[(7/5)(ln n+9)+10/9]<2(ln n+10)/n.

For j1 the sharper inverse-alpha estimate gives
|I1|<(25/18)(ln n+8)/n<2(ln n+9)/n.
The whole T(b)=-integral dv I1(A_b) therefore obeys
|T|<4(ln n+9)/n. The real spacelike contour obeys the same bound.
Independent complex u-integrations check all three full moments.

## Both whole matter-graviton endpoints

The light-active S290 F1 numerator integrates radially to I1-2I2+I3;
its F2 numerator gives[(2v-1)^2 I3-I1]/2. The contour satisfies
|2v-1|^2=r[1+(1-r)^2/4]<=1. The heavy-active denominator obeys
n u+(1-u)^2-a u^2(1-v^2)/4>=1+(n-5)u>=1 for a<=16.
Its moments and all zero-transfer OS terms can be bounded on the real
path. Since(ln n+10)/n is decreasing and ln128<5, it is<=15/128.
Writing c=g^2/(16pi^2), the whole light F1 plus the heavy and two
subtraction terms has |f1(s)|<c[16*15/128+1/4]<3c.
Thus |f1(s)/s|<12c/25<c. For a<=0 the exact OS difference quotient
has integrand weight*Aweight/(Delta_a Delta_0); both denominators
are>=1, so |f1(a)/a|<=c/2, continuously at a0.

The light F2 bound, heavy F2 and OS constant give
|f2_tri|<c[4*15/128+1/4+1/12]<c. The H-metric bubble adds
c|Bbar|/(n-a)<16c/n<=c/8, so |f2_g|<2c. No direct C bubble is
added here: that entire endpoint is already included in S293.
