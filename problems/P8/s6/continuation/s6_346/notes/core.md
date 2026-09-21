# Complete selected sum and exact tuned cancellation

Write c=C/g^2, a=c+1/n. The bubble expansion follows directly from

    B_MS(v)=v/6+v^2/60+v^3/420+O(v^4).

The v^3 coefficient of A(v)^2 B(v)/(2g^4) is

    L=a^2/840+a/(60n^2)+a/(6n^3)+1/(12n^4).

Consequently, before common g^4/(16pi^2),

    c_bubble=-4L,
    c_triangle=a J0+J1/n^2+J2/n^3+J3/n^4,
    c_core=c_bubble+c_triangle+c_box.

All four heavy-resolvent orders and all24 triangle labelings are included.
Their normalizations are inherited from the full exact covariant
functionals, not independently chosen effective vertices.

Substitute the ORIGINAL c=-3/(n-2)+2/(n-2)^2 exactly.
The entire rational-log expression is stored in the report.
Direct symbolic limits give n^2*c_core ->0, n^3*c_core ->0 and
n^4*c_core ->-23/105. These cancellations are a property of the
common-basis aggregate; separate sector bounds do not establish them.

The real-TT difference at this homogeneous order is
chi_core*T/sqrt(kappa), chi_core=g^4*c_core/(16pi^2).
The constant OS4 subtraction contributes no degree6 term.
Full physical amplitudes already contain these loops; adding this
contact on top of the unexpanded loops would double count.
