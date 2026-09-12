# The actual formal adjoint and its distinct causal inverse

The detector-source pairing uses conformal Lebesgue time density. For Bc=L D^2+Q+C D its formal adjoint is

Bc*=L^T D^2+Q^T-D(C^T)
    =B0^T-C^T D-(C')^T.

In particular, transposing the two channel indices is not enough; the derivative of the time-dependent coefficient is mandatory. The complete pointwise integration-by-parts boundary is

z^T L y' - z'^T L y + z^T C y,

whose derivative is z^T Bc y-(Bc* z)^T y. Production verifies this symbolic identity. Three independent compact-polynomial fixtures atq=0,7/5,10000 verify the integrated pairing with nonconstant h=1+t+t^2/3 and detect omission ofC'.

Let R0^T denote the transpose of the flat constant-coefficient causal inverse kernel. Seek Bc* z=f by

z=R0^T*f+R0^T*(C^T z')+R0^T*((C')^T z).

With v=z' and z=Iv, R0(0)=0 gives

v=(R0')^T*f+(R0')^T*(C^T v)+(R0')^T*((C')^T Iv).

Reorder the last double integral so that it is an integral againstv. Its effective kernel at(t,s) is the sum of(R0')^T(t-s)C(s)^T and

integral_s^t (R0')^T(t-u)(C'(u))^T du.

Thus its norm is at mostC0(8+35(t-s))<=C0(8+35T)<=215/2<108 forT<=1. Volterra's factorial estimate, applied separately to this actual adjoint equation, gives

||partial_t Z(t,s)||<=C0 exp(108(t-s)),
||Z(t,s)||<=C0(t-s)exp(108(t-s)),
Z(s,s)=0,

where Z=(Bc*)^-1 is the RETARDED inverse of the formal adjoint. Equation and uniqueness prove both inverse identities on the complete causal graph, as forY.

Z is not obtained by transposing the forward kernel Y(t,s) at identical ordered times: the actual operator coefficients and derivative order differ. Also, the alternative equation z=R0^T*f+(R0')^T*(C^T z) bounds the undifferentiated response but does not alone imply its derivative bound. The separate velocity proof above avoids that gap.

Finite lower-triangular time/channel matrix fixtures test both inverse products, deletion of the coefficient commutator and reversed factor order. They are explicit noncommuting algebra diagnostics, not certificates that a finite-difference discretization converges to the continuum operator.
