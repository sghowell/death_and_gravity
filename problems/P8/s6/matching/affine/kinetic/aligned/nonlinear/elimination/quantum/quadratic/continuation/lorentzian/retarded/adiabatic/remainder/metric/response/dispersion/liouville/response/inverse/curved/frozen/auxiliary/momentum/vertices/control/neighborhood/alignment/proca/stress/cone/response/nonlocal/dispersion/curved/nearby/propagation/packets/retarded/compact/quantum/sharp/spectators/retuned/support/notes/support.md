# Entire spatial multiplier and the actual matter cone

For zeta in C^3 use q=(zeta dot zeta)/R^2, a bilinear square, not a
Hermitian norm. On every compact zeta set the original polynomial
generator is uniformly bounded in time. Its Peano--Baker series converges
uniformly there, as do its parameter derivatives. Hence U_X and the
original local scalar response G_hat are entire functions of zeta.
This includes zeta dot zeta=0; no globally holomorphic choice of k is
required. Choose either square root k pointwise for the estimates only.

Write zeta=x+i*y. Exact real algebra gives

    (|x|^2+|y|^2)^2-|zeta dot zeta|^2
      =4*sum_{i<j}(x_i*y_j-x_j*y_i)^2 >=0.

Consequently |k|<=|zeta| and
|Im(k)|^2=(|k|^2-|x|^2+|y|^2)/2<=|y|^2. The large complex null vector
A*(1,i,0) illustrates why a small scalar-root disc is not a small
complex-vector ball. It is explicitly included by the original chart.

Combining the COMPLETE high- and low-root estimates, uniformly for
s<=t in I, gives

    |G_hat(t,s,zeta)| <= (4+320*exp(400))*exp(S_m(t,s)*|Im(zeta)|).

The radius here is the actual integral S_m, not its convenient numerical
upper estimate. The identity -N^2+(eR)^2*(N/(eR))^2=0 checks that this
integrand is exactly the physical radial null coordinate speed.

The distributional Paley--Wiener--Schwartz theorem turns this entire
exponential-type bound into support in the closed radius-S_m ball.
For the sharp radius, pair with a test beyond a separating hyperplane
and shift the Fourier contour toward its normal. The multiplier's
growth is exp(eta*S_m), while the test contributes
exp(-eta*(S_m+d)) times polynomial factors, d>0; the pairing vanishes.
Rotating the normal gives the ball. See the
[MIT distributional theorem and contour exercises](https://math.mit.edu/~rbm/Problems4.pdf).

The real-momentum multiplier is uniformly bounded and continuous in
(s,t), so its inverse Fourier transform is a continuous tempered
distribution-valued function on ordered times. Multiplication by the
retarded time ordering defines the spacetime response. S_m is continuous
and additive along this FLRW solution; its balls are exactly the
physical-matter causal reach. Compact spacelike separated source and
detector tests consequently have zero COMPLETE retarded pairing.

At equal times U_X is the identity and its chi-to-P_chi entry is zero.
No instantaneous tail is omitted. This proves support, not a uniform
interacting error bound or a statement about a different parent theory.
