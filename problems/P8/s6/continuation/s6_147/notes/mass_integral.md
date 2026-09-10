# Regulator-first outer mass integral

For general 0<b<=1 the outer triangle is

    T3_D(v,b,b)=exp(gamma_E e)mu^(2e)Gamma(e)/Q
       [b^(-e)/(v-b)
        -(v^(1-e)-b^(1-e))/((1-e)(v-b)^2)].

A Feynman-parameter integral or the radial integral
Integral_0^infinity x^(1-e)/[(x+v)(x+b)^2] dx / Gamma(2-e)
gives this formula. Its four-dimensional limit is
[v log(v/b)-v+b]/[Q(v-b)^2]. The apparent v=b singularity
is removable. The actual domain v>=4m^2 is far from it.

Use the complete regulated signed density and proper MS subtraction
of S6.144, not a positive spectral-measure assumption. Let ell=log(m^2).
Keeping both leading soft and hard terms, with CY=2NY^2/Q, gives

    F_mix,D(b)/(CY/Q)
      = exp[(ell-log b)e]
          [(12-32e+16e^2)A(e)^2-12A(e)]/e^2
        -6 B(e)/e^2 + delta_D(b),

where

    B(e)=sqrt(pi) exp(2 gamma_E e)4^(-e)(1-2e/3)(1-4e)
         Gamma(1+e)Gamma(1+2e)
         /[(1-e)(1+2e)Gamma(1/2+e)].

This hard term is independent of b. The finite leading sum is
-46-32(ell-log b). The full difference delta_D(b) is convergent
before e=0 and is treated separately in remainder.md.

Since half of CY/Q equals NY^2/Q^2, the normalized mass difference
is minus the b integral. Integrating at nonzero regulator gives
exactly

    {-exp(ell e)[(12-32e+16e^2)A(e)^2-12A(e)]/(1-e)
       +6B(e)}/e^2
       =6/e^2-2/e+78+32ell+O(e).

This is independently equal to the integral of the finite leading
coefficient, because Integral_0^1(-log b) db=1. The agreement is
between complete regulated forests; it does not license an early
d=4 subtraction of either raw piece.

For |e|<1/8, multiply by e^2 to remove the displayed poles.
The remaining leading functions have an integrable b^(-1/8)
majorant with any fixed finite power of |log b| for their e jets.
The difference integral has an additional b/T suppression,
T=4m^2. The signed-density endpoint powers are integrable on this
complex neighborhood and the ultraviolet difference has an extra
1/v. Uniform domination permits coefficientwise b integration.
The overall MS operation then removes the two displayed poles.

An independent low-energy check is the quartic threshold
Gamma_0,MS(0)=-64NY^2/Q. Half its light tadpole is
32NY^2(ell+1)/Q^2, reproducing the log coefficient 32.
The remaining leading constant 46 is a hard-region contribution,
not an arbitrary new local parameter.
