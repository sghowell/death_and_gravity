# Actual composite coordinate: bounded window, nonuniform extension

Use the exact S6.118 MSbar composite prescription, with c=4lambda^2/gamma.
At Euclidean Box eigenvalue y>=0 the finite mixing is N(y)/(16pi^2),
where

    N(y)=c/2-3lambda+gamma/2
          -(lambda+5gamma/4)y + gamma y^2/2.

Native code rebuilds this from the frozen actual polynomial.
Its completed square is

    N(y)=gamma/2 [y-lambda/gamma-5/4]^2
          +3lambda^2/(2gamma)-17lambda/4-9gamma/32.

The minimum term is strictly positive at the actual parameters.
On 0<=y<=M, coefficientwise absolute values and pi>3 give

    |W_MS(y)| < [c/2+3lambda+gamma/2
                 +(lambda+5gamma/4)M+gamma M^2/2]/144
              < 10^-404.

This window includes every real Euclidean momentum up to the
heavy mass. It is not a restriction imposed on loop integration.

The same polynomial grows quadratically at infinity. An explicit
finite witness is y=10^400: it lies above M but below the selected
Planck mass squared 10^800. Exact arithmetic gives N(y)>256;
pi<4 then gives W_MS(y)>1. No floating cancellation is used.

Hence this specific one-loop composite coordinate does not have
a uniform small relative kinetic correction on the full massive
H1 domain. Fourier wave packets normalized in that domain and
localized at increasingly large momenta exhibit the unbounded
multiplier directly.

This is a statement about the specified off-shell composite
normalization. Finite composite counterterms can change that
off-shell polynomial. The witness is not a physical cutoff,
a ghost, a failure of stable-light LSZ equivalence or an
obstruction to another field representation. In particular it
does not invalidate S6.118's proper on-shell amplitude transfer.
It prevents using the low-momentum norm estimate as an
unqualified global off-shell estimate.
