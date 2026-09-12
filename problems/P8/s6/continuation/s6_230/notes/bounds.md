# Explicit massive-factor bounds and actual pole scale

All logarithms are on the ORIGINAL first sheet. Write r=|p|/m^2, delta=4m^2/p, d=1+delta, s=sqrt(d), P(d)=30-20d+3d^2. For r>=8, |delta|<=1/2. The original exterior identity from S226 is

atanh(1/s)=[Log(p/(4m^2))+2Log(1+s)]/2.

Thus, without replacing the massive form factor by a massless approximation,

A2=B(d)+T(d)[Log(p/m^2)+Lc(d)],
B=-172/225+19d/30-d^2/10,
T=s P(d)/60, Lc=2Log[(1+s)/2].

These branches agree on positive p and extend analytically throughout the exterior first sheet, including their limiting cut banks. Set alpha=13/60, beta=-52/225.

## Uniform error, including near the cut

The square root has positive real part. On |delta|<=1/2,

|s-1|<=|delta|, |s|<5/4, |1/s|<3/2, |1+s|>1.

Also P(1+delta)-13=delta(-14+3delta), so its modulus is at most16|delta| and |P|<=21. Consequently

|T|<1/2, |T-alpha|<=37|delta|/60<|delta|,
B-beta=delta(13/30-delta/10), |B-beta|<|delta|/2.

Since (1+s)/2 differs from1 by at most |delta|/2<=1/4, the convergent logarithm series gives |Lc|<=4|delta|/3. Combining these bounds and |Log(p/m^2)|<=log r+pi yields

|A2-alpha Log(p/m^2)-beta|<=8(1+log r+pi)/r.                    (1)

The same estimates imply the convenient loose exterior bound |A2|<=20+2 log r. On the outer circle r8 this is below26. On the finite cut segment tau/m^2 in[4,8], put v=1-4m^2/tau in[0,1/2]. The real local part of D_A lies between-172/225 and0, while the positive atanh term is at most1, using atanh(sqrt(v))<=sqrt(v)/(1-v). Hence |D_A|<2. Also U<=1/2, so |A2|<4 on either bank, using pi<4. A2 is finite at threshold and removable at p0. The maximum principle on the slit disc, followed by its boundary limit, gives |A2|<30 for r<=8. Together,

|A2(p)|<=30[1+log(1+|p|/m^2)]                                (2)

on the complete first sheet. The branch point is handled by a vanishing indentation; no singular boundary is skipped.

## Certified coarse enclosure for the actual C

Here m=1000, kappa=10^800, C0=16pi^2 kappa and C=C0+5m^2/6. Elementary 3<pi and pi^2<10 give

144 kappa<C<161 kappa.

At a complete-denominator root z, |z A2(z)|=C. If r=|z|/m^2<=10^790, (2) gives |A2|<30(2+790*3)<10^5, whereas C/m^2>10^796. This is impossible. Therefore |z|>10^796.

For this exterior root, the original S226 minorant is

Re A2>=13 log(r/4)/80-4>250,

using log10>2 and log4<2. It follows that |z|=C/|A2|<C/250<kappa, and hence r<10^794. Thus the rigorous enclosure is

10^796<|z|<10^800,
10^398<|sqrt(z)|<10^400.                                     (3)

This is between one hundredth and one times sqrt(kappa), not ABOVE that normalization.

For the upper root, theta=arg z is in(pi/2,pi). At its radius, (1) has error below10^-785, since 8(1+794*3+4)<10^5. Therefore

Im A2>alpha*3/2-10^-785>3/10, |A2|<20+2*794*3<5000.

Since z=-C/A2, Im z=C Im A2/|A2|^2>kappa/10^6. The positive-real square root obeys

Re sqrt(z)=Im z/sqrt[2(|z|-Re z)]>sqrt(kappa)/(2*10^6)>10^393. (4)

Every rational margin in these coarse comparisons is replayed. Numerical pole approximations in the private diagnostics are not substituted for (3)-(4).

## Residue bound without crossing a cut

Differentiate the original analytic coefficient formula directly. Because p delta'=-delta,

p A2'=T-delta[B_d+T_d(Log(p/m^2)+Lc)+T Lc_d].

Here |B_d|<1; the derivative formula for T gives
|T_d|< [21*(3/4)+(5/4)*17]/60<1; and
|Lc_d|=|1/[s(1+s)]|<3/2.

At the actual root, the above radius/logarithm bounds make the delta correction smaller than1/2, while |T|<1/2. Hence |z A2'(z)|<1. No Cauchy disc crossing the cut is used. With Re A2>250,

R=1/[z D_C'(z)]=-1/{C[1+z A2'(z)/A2(z)]},
|R|<=250/(249 C)<2/C.                                      (5)

The strict bounds and simplicity from the complete count both exclude a zero residue.

## Original low-energy physical normalization

For |p|<=m^2 the defining denominator has modulus at least3m^2 and integral W2=3/14. Thus

|A2|<=1/30+(3/14)/3=11/105.

There are no additional full-response poles in this disc; the Einstein massless p0 pole remains. Compare the physical propagator to the Einstein propagator through their ratio, including the continuous ratio at p0:

|C0/(C+pA2)-1|
 <=(5m^2/6+11m^2/105)/(C-11m^2/105)
 <197m^2/(30240 kappa)<10^-796.                              (6)

The original finite m^2 R shift is not dropped or retuned. This is a complex Lorentz-invariant form-factor estimate, not a temporal bandpass theorem, a physical EFT cutoff, or protection against exciting the distant poles by compact-time forcing.
