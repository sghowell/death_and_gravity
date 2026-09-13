# Complete four-order loop and full logarithmic moment reduction

The exact S236 primitive is T=w[P(w)H0+Q(w)], with w=1/n and H0=ln n-Log L. The complete analytic coefficient jets through degree3 are

P0=1, P1=3L, P2=10L²-2Lb, P3=35L³-15L²b,

Q0=-1,
Q1=(b-7L)/2,
Q2=(-74L²+28Lb+b²)/6,
Q3=(-533L³+327L²b-9Lb²+b³)/12.

The logarithm appearing in Q is log((1+d)/2)=log(1+(d-1)/2), where d(0)=1. Its log1p increment is(d-1)/2, not(1+d)/2. Independent derivative checks retain this distinction.

For L=1-s xi(1-xi), let Mj=integral L^j and Jj=integral L^j Log L. The complete first moments are

M0=1,
M1=1-s/6,
M2=1-s/3+s²/30,
M3=1-s/2+s²/10-s³/140.

With y=xi-1/2 and a0=1-s/4, integrate the FULL derivative of y L^j Log L. Its endpoint terms vanish because L=1 at y=+/-1/2. For j>=1,

(2j+1)Jj-2j a0 J(j-1)+2Mj-2a0 M(j-1)=0.

This is valid throughout the specified complex bidisk, where Re L>0 and no logarithmic branch is crossed. It reduces the logarithms to J0 and explicit polynomials; the Jj terms are not set to zero.

For the second HEAVY box parameter b=t eta(1-eta),

integral b^k=t^k(k!)²/(2k+1)!.

The triangle averages[Pj H0+Qj]/n^(j+1) at b0. The complete box is its full negative heavy-mass derivative and averages

[(j+1)Pj H0+(j+1)Qj-Pj]/n^(j+2).

The last -Pj is from differentiating ln n and is essential. The exact vertex identity is

A_s/g²=-2w+(s-4)w²+(s²-4)w³+s³w4+r_a,
r_a=w5[s4/(1-sw)+16/(1-2w)²].

Combine every channel bubble, all triangle terms and all six ordered boxes. Only after imposing s+t+u=4 and the full logarithmic moment identities does the J0 coefficient vanish in each light channel through n^-5.

Writing sigma2=s²+t²+u² and sigma3=stu, the full amplitude before g4/(16pi²) is

F2/n²+F3/n³+F4/n4+F5/n5+remainder,

F2=-6ln n,
F3=-16ln n+10/3,
F4=(3ln n/2-157/90)sigma2-28ln n+10/9,
F5=(3ln n-13/21)sigma3+(12ln n-7969/630)sigma2-28ln n-2951/63.

The first two orders are constant; the next orders have total degrees2 and3. A purely formal series would not bound coefficient errors. The following notes give the complete analytic remainder instead.
