# Full finite-width derivatives, not a clock substitution

Use exactly the S263 real lapse strip[999/1000,1001/1000]
and its entire-source C5 comparison envelope below10^-380.
The original clock comparison is R_clock=1+(N^-2-1)/h.
Its first four lapse derivatives have absolute bounds3,7,25,121.
Adding the FULL source error gives safe original derivative bounds
4,8,26,122. The full local R also lies in[1/2,2].
Smooth fixed time profiles remain those already controlled; no
complex-analytic assumption about them is added.

For each alpha=1/4,1/2,3/4 and derivative order0..4, factor
d_N^j(R^-alpha)=R^-alpha P_j(R,R_N,...,R_NNNN).
Each P_j is a rational Laurent polynomial. Exact symbolic
reconstruction and outward rational interval arithmetic give
the derivative bounds
[2,4,48,986,29002],
[2,8,112,2548,80692] and
[2,12,192,4782,161982], respectively.
All are below the safe common ceiling10^8.

The complete original difference from the clock has order8 and
order1024 source factors at X=1. It therefore licenses the specified
finite CENTER jets, not replacement of the whole off-center function.
The exact first three derivatives are
c1=2alpha/h,
c2=4alpha(alpha+1)/h^2-6alpha/h,
c3=8alpha(alpha+1)(alpha+2)/h^3
 -36alpha(alpha+1)/h^2+24alpha/h.
All five volume-factor jets are checked against the original S253
lapse-jet source, preserving the S262 physical binding.

At the bounce, alpha=3/4 gives jets
1,3/2,3/4,-3/8 through third order.
The FULL function is used outside the Taylor event |n|<=10^-3.
The Gaussian tail is not set to zero, and local derivative bounds
are not assumed valid on that unbounded complementary event.
