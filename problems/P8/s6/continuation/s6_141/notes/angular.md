# Exact four-dimensional angular integral

Write x=q^2 and y=l^2. Rotational invariance
of the positive norm majorant leaves the
normalized relative-angle measure
(2/pi) sin(theta)^2 dtheta on [0,pi].
For a=x+y, b=2sqrt(xy), a>b>0, polynomial
division gives

 (1-c^2)/(a-bc)
   =c/b+a/b^2+(1-a^2/b^2)/(a-bc).

The first term integrates to zero. In the
last term use t=tan(theta/2):

 integral_0^pi dtheta/(a-b cos theta)
   =2 integral_0^infinity dt/[(a-b)+(a+b)t^2]
   =pi/sqrt(a^2-b^2).

Thus the normalized chord angular average is

 2[a-sqrt(a^2-b^2)]/b^2 =1/max(x,y).

The square root is |x-y|, so both radial orders
give the same piecewise formula. When x=y>0,
the numerator sin(theta)^2 cancels the
small-angle divergence of the denominator;
the integral equals 1/x directly and by its
continuous limit. When one radius vanishes the
denominator is the other radius squared. The
single point x=y=0 is handled by the subsequent
finite positive double radial integral.

Adding the scalar mass b_scalar=1 can only
decrease the positive unshifted chord. The
same massless expression bounds both sectors.
No angular average of a complex-shifted or
sign-indefinite propagator is asserted.
