# The actual source and sufficient light-field derivative bounds

The S6.42 expansion retained Q2 symbolically. It is not a new free
coupling. Differentiate the original Q_x+A Q=F equation at x=-1:
Q=Q_x=0 implies Q_xx=F_x=-3h_phi/(4h³). Since x+1=2eps*n+O(eps²),

    Q2=2Q_xx*n²=-9H*n²/(4h²),
    S2=A*n*delta_K+C*n²,
    A=-2/h, C=-6H/h-27H/(8h²).

This substitutes into the exact full stationary source, not a newly
postulated driving term. In unitary gauge phi is unperturbed and the
background coefficients depend only on time. Metric and vector mass
coefficient changes multiplying a second-order heavy response are of
higher order and are not included in this leading estimate.

## Compact coefficient bounds

On |u|<=1/2, exact numerator coefficient sums and positive even
denominators bound the four coefficient derivatives as follows:

| Derivative order | Bound for A | Bound for C |
|---|---:|---:|
| 0 | 2 | 483/16 |
| 1 | 6 | 4497/32 |
| 2 | 33 | 3729/4 |
| 3 | 228 | 8982 |

For each differentiated rational function, leading.py verifies that
its denominator has positive constant and nonnegative even coefficients.
Its numerator is bounded by sum_j |coefficient_j|/2^j and its denominator
by its constant term. Thus these are continuous compact-domain bounds.

If all time derivatives through order three of n and delta_K are bounded
by epsilon_n and epsilon_K, Leibniz gives source-jet bounds

    mixed:        [2, 10, 65, 514]*epsilon_n*epsilon_K,
    squared lapse:[483/16, 6429/32, 12921/8, 132027/8]*epsilon_n².

The common sufficient bound is

    E=514*epsilon_n*epsilon_K+(132027/8)*epsilon_n².

It is a bound for the actual quadratic source and its first three time
derivatives, not a claim about the full nonlinear source. Preparing
n=n'=n''=0 at the left endpoint gives S2=S2'=S2''=0 there without
requiring delta_K to vanish. The on-shell light evolution has not been
proved to maintain these derivative bounds or to realize every profile.

## Products are Fourier convolutions

For the mode statement, use smooth periodic perturbation profiles as
Fourier bookkeeping on the unchanged comoving physical background.
No change of spatial topology or action is needed. Let

    e_n(k)=max_(j=0..3) sup_interval |partial_u^j n_k|,
    e_K(k)=max_(j=0..3) sup_interval |partial_u^j delta_K_k|,
    sum_k e_n(k)<=epsilon_n, sum_k e_K(k)<=epsilon_K.

Absolute summability legitimizes multiplication and the stated time
derivatives. For each output mode m, all four source jets are bounded by

    E_m=514*(e_n*e_K)_m+(132027/8)*(e_n*e_n)_m,
    sum_m E_m<=E.

The stars are convolution, not pointwise multiplication of coefficients.
In particular, a single incoming cosine sources zero and double momentum.
The response theorem applies to each output mode with 0<|m|²<=1,
using its E_m; it gives the summed physical low-mode bounds with E
as well. The homogeneous source is an algebraic temporal control.
No estimate is made for output momenta outside this declared window.

Complex Fourier coefficients obey the same energy proof with absolute
squares and real parts of inner products. A real field supplies the
usual conjugate pairs. The vector spatial readout uses the norm of
its unit momentum direction, so summing mode bounds uses the ordinary
triangle inequality with no hidden polarization factor.

The time-jet and summability assumptions are stronger than the earlier
spatial source bound. They are explicit sufficient assumptions, not a
derived nonlinear evolution theorem, a physical frequency band, a
Wilsonian cutoff, or an original P8 completion condition.
