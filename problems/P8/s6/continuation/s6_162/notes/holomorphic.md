# One common complex field-amplitude circle

Let all 21 base canonical jet components have modulus at most
one, and scale all of them by a common complex t. This is field
amplitude scaling at fixed jets, not a momentum scaling.
Choose rho=10^300. On |t|<=rho,

    |u|<=U=10^-100, |X|<=Xcap=4*10^-200.

The following bounds apply to the full rational T and full R.
For |X|<=1/4 let r=X/(1-X), so |r|<=1/3 and
T=r^n/(1+r^n). Thus

    |T|<=2*3^-n,
    |T/X|<=8*3^-n,
    |T'/X|<=100n*3^-n.

For example T'/X=n r^(n-2)/(1-X)^3/(1+r^n)^2,
whose norm is strictly below the displayed cap.
At X=0 use the removable analytic values. The exact integer
inequality 2*10^450<3^1024 makes the first bound below 1e-450.

Since n Xcap^2<1/2, the exponential norm is below two and

    |w|<=2n Xcap^2,
    |w/X|<=2n Xcap,
    |w'/X|<=8n.

Writing these caps as Tup,Wup,TXup,WXup,TpXup,WpXup,

    Bup=Tup+(1+Tup)Wup,
    Rdef=(1+Xcap)Bup/(1-U^2)^3 <1/2,
    BXup=TXup+(1+Tup)WXup,
    BpXup=TpXup(1+Wup)+(1+Tup)WpXup,
    A3up=[(1+Xcap)BpXup+BXup]/(1-U^2)^3.

These give |R-1|<=Rdef, a nonzero R and

    |A3|<1e6,
    |A4|<=A3up+7(Xcap A3up)^2/[4(1-Rdef)]<2e6,
    |A5|<=A3up^2 Xcap/(1-Rdef)<8e-188.

The literal retuned tree numerator's coefficient L1 sum,
divided by 200(15/16)^12, is exactly

    354778132704900677632/3243658447265625 <1e8.

This bounds the tree throughout |u|,|X|<=1/4.
The actual Fv modulus is bounded by

    (Xcap+U^2)/2+(lambda kappa+n|F0|)Xcap^2
      +(n/3+|F0|)U^4 <1.

With |exp(-u^4)|<3 and |S|<4 this bounds |F| below 1e10.
Use |Y|<=4rho^2, |L3|,|L4|<=64rho^4, |L5|<=256rho^6.
All scalar and DHOST groups together have norm strictly
below M=10^811 for the canonical density on this circle.

All denominator and norm gaps are strict. The functions are
therefore holomorphic on a neighborhood of the closed circle,
uniformly for the entire base jet cube. The tiny but nonzero
rational step is never dropped in this full-function bound.
