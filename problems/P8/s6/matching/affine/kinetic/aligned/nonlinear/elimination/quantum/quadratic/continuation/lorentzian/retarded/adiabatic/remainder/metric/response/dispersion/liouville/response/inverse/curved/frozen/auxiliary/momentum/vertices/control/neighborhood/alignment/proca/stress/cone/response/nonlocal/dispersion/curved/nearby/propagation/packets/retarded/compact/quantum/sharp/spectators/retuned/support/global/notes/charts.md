# Complete overlapping Euler charts

The unitary chart uses x=(v,chi) and n=(v'+ell*chi/2)/Theta. Substitute
this into the complete literal light Lagrangian. Write the result as
a^3*(x'^T A x'/2+x'^T B x+x^T C x/2). The matrix A is q-independent;
B and C have degree at most one in q, and the q coefficient of B is
symmetric. The exact equation is

    x''+L*x'+Q*x=0,
    L=A^-1*(A'+3H*A+B-B^T),
    Q=A^-1*(B'+3H*B-C), q'=-2Hq.

The gamma chart uses x=(b,chi) and the pinned moving canonical exchange
P_b=2*a^3*q*v. The Hamiltonian retains -H*b*P_b. Its EXACT Legendre
matrices A,B,C are rational in q; they are not truncated to their
principal limits before the Euler derivative. Both charts satisfy

    K=lim A/2, G=lim(B'+3H*B-C)/(2q),
    K=K0+diag(delta_J/d^2,0),
    K0=[(J+w^2/2)/d^2, sign*w/(2d); sign*w/(2d),1/2],

with (d,sign)=(Theta,+1) or (Lambda,-1). The ACTUAL identities
w=-ell*Lambda and
J+w^2/2=Theta*(H*Lambda+Lambda')-Lambda*Theta'-Theta^2
give G=K0 in both charts. These identities are replayed from the actual
global background, not imposed on independent arbitrary coefficients.

The complete friction L and remainder Q-q*K^-1*G have rational q-degree
at most zero. They are exactly q-independent in the unitary chart.
Native coefficientwise degree checks and exact full 4x4 intertwining
identities verify the Euler equations against the regular density
generator, including every time jet and the moving boundary. The maps
from Z to Y=(k*x,x') have frequency degree at most two, and their exact
inverses have degree at most one.
Every denominator factor of those maps and full remainders is checked
against the declared pivots: a,k,J_new and Theta in the unitary chart;
a,k,J_new,Lambda and the gamma D factor in the gamma chart.

In the gamma chart the finite-q denominator is
D=q*Lambda^2-J-delta_J-w^2/2. At the bounce it vanishes at q=152/25.
The original generator is regular there. We use the gamma velocity
description only at the explicitly separated large COMPLEX frequencies;
lower frequencies are always evolved in the original phase.

The all-time chart cover is explicit. On |u|<=1/4,
h<5/4, -1/2<=Lambda<-1/5 and a^2<(17/16)^4<2. Since |ell|<=1/10,
w^2/2<=1/800; the earlier |u|<=1/2 bound gives J<36. Thus
J+delta_J+w^2/2<36+1/50+1/800. If |k|>=64, then |q|>2048 and
|D|>=|q|*Lambda^2/2>0. This separates every complex gamma pole.

On |u|>=1/8 within a finite strip |u|<=R, R>=1/4,
|Theta|>=3/[8*(1+R^2)]>0. The unitary chart has no frequency pole.
Switching at u=+-3/16 stays inside both overlaps and divides ANY finite
ordered interval into at most three regular segments. No crossing is
skipped, no coefficient is reset and no global uniform-in-time norm
constant is presumed.
