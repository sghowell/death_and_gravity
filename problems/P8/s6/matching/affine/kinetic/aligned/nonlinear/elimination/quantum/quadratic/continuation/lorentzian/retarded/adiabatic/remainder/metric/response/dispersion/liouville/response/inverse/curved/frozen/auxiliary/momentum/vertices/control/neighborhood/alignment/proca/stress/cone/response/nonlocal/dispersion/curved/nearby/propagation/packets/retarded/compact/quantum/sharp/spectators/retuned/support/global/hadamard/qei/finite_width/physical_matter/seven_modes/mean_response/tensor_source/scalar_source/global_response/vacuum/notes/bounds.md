# Exact whole-tube mixed norms and uniform order control

On 9/10<=X<=11/10, all real u, put delta=21/100.
Use the finite normalized mixed norm

    ||f||_4=sum_(a+b<=4) sup |d_u^a d_X^b f|/(a! b!).

For curvature/DHOST differences apply the external weight h
to each derivative; for the scalar difference apply h^2.
These are derivatives of f followed by a weight, not
derivatives of h f or h^2 f. Weighted product estimates follow
by Leibniz with normalized factorials.

For j<=6, the exact quadratic-base differentiation recurrence
gives the majorant

    C_j=sum_(b=0..floor(j/2))
      j! (n)_(j-b) delta^(n-j+b) (11/5)^(j-2b)
      /[(j-2b)! b!].

The falling factorial and all terms are positive. The native
checks independently verify the recurrence polynomials.
Let S=sum_(j<=4) C_j/j!, Sd=sum_(j<=4) C_(j+1)/j!.
The normalized inverse-h derivative sum is at most
H=sum_(a<=4) (6)_a/a!=210. To prove the weighted version,
factor h^-1=(u-i)^-3(u+i)^-3. Each extra derivative produces
a denominator factor of magnitude at least one, while h
cancels the base product. Vandermonde gives (6)_a.

Take Xm=11/10 for ||X-1|| and Xinv=sum_(j<=4)(10/9)^(j+1).
Write E=S_n+(X-1)S_n' and F=(X-1)S_n. Their norm bounds
are S+Xm Sd and Xm S. Rold=1+h^-1(X-1), Rnew=Rold-h^-1 F.
The native bound gives ||Rold||<=232, ||Rnew||<233 and
Rold,Rnew>4/5 on the closed clock tube. The mixed inverse
jet norm is bounded by

    (5/4) sum_(j=0..4) [(5/4)*233]^j.

This is the finite Taylor-algebra inverse at each point:
only four derivative orders contribute. It does not require
convergence of an infinite Neumann series.

With T=(-2E+E^2)Rold+h^-1 F, the exact differences are

    Delta A3=-h^-1 E/X,
    Delta A4=-Delta A3-(7/4)h^-2 T/(Rold Rnew),
    Delta A5=h^-2 T/(X Rold Rnew).

The displayed module retains every product in these bounds.
All h-weighted mixed C4 errors for F2,A3,A4,A5 are strictly
below 10^-400 at n=1024, kappa=10^800.

For the scalar, strip exp(-u^4) from successive derivatives:
P_(a+1)=P_a'-4u^3 P_a. Each rational denominator has positive
constant and nonnegative even coefficients, or the envelope
rejects it. Bound every numerator monomial using
|X|<=11/10 and |u|^p exp(-u^4)<=ceil(p/4)!.
Multiplying the resulting h^2-weighted mixed norm by S
bounds the full scalar difference, again below 10^-400.
This includes the calibrated u^4 potential and every
derivative of the field localizer. No time grid is used.

For the entire even-n co-scaled family, split exactly

    Fv_n-Ftree=n L+C,
    L=(lambda/gamma-F0)X^2-u^4/3,
    C=X/2-u^2/2-F0 u^4-Ftree.

Use the positive majorant n||exp(-u^4)L||+||exp(-u^4)C||,
rather than presuming cancellations survive varying n.
Its n=1024 bound still gives a scalar error below 10^-400.

On increasing n by two, each falling-factorial factor changes
by at most Lstep=1+2/(1024-5), and delta^n contributes delta^2.
All switch jets through order six contract by at most
delta^2 Lstep^6. The affine scalar envelope adds at most one
further Lstep. The resulting ratio is less than 1/20.
Positive sums and products in the curvature majorants also
contract; the inverse-R majorant 233 remains available.
Thus the stated four-jet budget holds for every even n>=1024.
It is a classical coefficient budget, not a loop-error bound.
