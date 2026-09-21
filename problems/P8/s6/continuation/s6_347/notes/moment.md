# Exact mixed bubble moments

For n>0, M=1+(n-1)z, set
b_r=(1/r)int_0^1[z(1-z)/M]^r dz for r2,3.
Substitute y=M and integrate every numerator monomial over y^r,
using log(y) for exponent-1. Check the full primitive derivative,
both endpoints and the removable equal-mass limit.

    b2=(n^3-6n^2logn+9n^2-6nlogn-9n-1)/(6(n-1)^5),
    b3=(n^4-12n^3logn+28n^3-36n^2logn-12nlogn-28n-1)/
        (12(n-1)^7).

The equal-mass values are1/60 and1/420. Combine all source classes:

    d=(n-1)b3-b2, E(n)=-32[d+2/n^4],
    chi_extra=g^2 E(n)/(16pi^2 kappa).

The exact expression and n^2E->8/3 are stored in the report.
At n1, E=-952/15; do not extend the original large-mass positivity
claim to arbitrary positive masses. Exact APIs retain this limit
and reject floating masses, nonpositive masses and invalid orders.
