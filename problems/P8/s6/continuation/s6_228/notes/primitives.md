# Ordered zero-past primitives and the weak logarithmic majorant

Translate the original proper-time interval to[0,1] only for causal integration notation. The background coefficient is still the original function evaluated at t-1/2; no physical time profile is changed. Every source and response uses the same nonempty zero initial neighborhood.

For j<=3, integrate the original local termc_j(t) h^(j)(t) four times in the OUTPUT time. Integration by parts on its actual source variable gives

I4(c_j h^(j))(t)=integral_0^t K_j(t,s)h(s)ds,
K_j=sum_(r=0)^j (-1)^r binom(j,r)
    (t-s)^(3-j+r)c_j^(r)(s)/(3-j+r)!.

All exponents are nonnegative. The upper boundary vanishes before the final ordinary integral because the original primitive has degree3 and only j<=3 source derivatives are moved. All lower boundaries vanish by the common initial germ. The coefficients are differentiated at the source time; they are not pulled through the primitives.

Independent exact polynomial integrals verify all four identities. For everyj>0, retaining only the undifferentiated coefficient changes that integral. These are checks of the ordered local formulas, not a numerical approximation to the full quantum kernel.

Together with the fixed singular extensions in notes/remainder.md, the actual total remainder has an ordinary kernel with bound

|V(t,s)|<=C(1-log(t-s)), 0<t-s<=1,

for some finite positiveC. It includes the entire actual quantum remainder, original contacts, known local lower coefficients and the classical terms proportional to the unchanged kappa.

For lambda>=1, substituting y=lambda r gives

integral_0^1 exp(-lambda r)(1-log r)dr
 <=(2+log lambda)/lambda
 <=2/sqrt(lambda).

For the first inequality, bound the constant part by1+log lambda, bound integral_0^1 -log y dy by1, and drop the nonpositive -log y contribution fory>=1. For the second, putu=sqrt(lambda)>=1 and use logu<=u-1. Equivalently, the margin2u-2-2logu vanishes atu1 and has derivative2(u-1)/u>=0.

Weighted Young's inequality therefore gives ||V||<=2C/sqrt(lambda) on continuous functions with norm sup exp(-lambda t)|h(t)|. This bound is on the actual constructed ordinary kernel. It is not inferred from the differential order of an unbounded operator or from the older weak response estimate.
