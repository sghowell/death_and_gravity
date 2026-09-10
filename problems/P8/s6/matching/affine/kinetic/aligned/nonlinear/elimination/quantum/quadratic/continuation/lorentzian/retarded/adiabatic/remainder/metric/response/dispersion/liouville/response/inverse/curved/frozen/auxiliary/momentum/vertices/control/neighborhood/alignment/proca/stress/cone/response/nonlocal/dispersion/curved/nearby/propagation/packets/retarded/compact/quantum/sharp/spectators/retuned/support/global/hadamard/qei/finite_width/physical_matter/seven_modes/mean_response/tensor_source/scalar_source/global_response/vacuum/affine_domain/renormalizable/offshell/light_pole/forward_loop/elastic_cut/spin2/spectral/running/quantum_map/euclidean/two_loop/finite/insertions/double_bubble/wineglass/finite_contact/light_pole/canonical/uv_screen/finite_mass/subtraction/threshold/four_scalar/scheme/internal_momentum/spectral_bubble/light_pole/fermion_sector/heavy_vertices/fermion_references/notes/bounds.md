# Rational reference enclosures

The moment bounds give J0+2R in [0,2]. With
upper bounds Yhi, ahi, Q>=Qlo>0 and Cf=4/3,

    0<=z<=zhi=(3Yhi/4+ahi Cf/2)/Qlo,
    |eta|<=eta_abs=(Yhi+2ahi Cf)/Qlo,
    |upsilon|<=ups_abs=(2Yhi+6ahi Cf)/Qlo.

Obtain w_abs from the SAME complete physical
Phi kinetic interval in S6.133; do not reset
that normalization at zero momentum.

For |w|<=1/2 and z>=0, let t=sqrt(1+w)>1/2.
Then 1/t<2 and
|1/t-1|=|w|/[t(t+1)]<=2|w|. Direct algebra
therefore gives, for the selected h=1 ratios,

    |m_zero/m_MS-1|<=eta_abs+zhi,
    |y_zero/y_MS-1|<=2(ups_abs+zhi+w_abs).

The corresponding first-order Yukawa bound
is ups_abs+zhi+w_abs/2. If eta_abs<1 and
ups_abs<1 the selected mass and Yukawa
reference ratios are positive. Large valid
inputs can instead return inconclusive
positivity flags; the API does not falsely
certify them.

For the actual candidate the positive quartic
ray lies between 3/2 and 2, giving ahi=2L/3.
The inherited rational Yhi and Qlo=144 yield

    zhi approximately 1.04197749385 times 10^-207,
    mass-ratio error <=4.31508193781 times 10^-207,
    Yukawa-ratio error <=2.79862203716 times 10^-206.

All three are below 10^-204 by exact rational
comparison. These are local-reference shift
bounds, not estimates of the remaining
two-loop primitive integrals or a fermion
pole spectrum.
