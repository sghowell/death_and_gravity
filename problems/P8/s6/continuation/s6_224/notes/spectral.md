# Shifted trace and shear measures with uniform primitives

The source-pinned scalar factors are Ftrace=-Atrace and F2=-A2. Their positive-sign radial weights are Wtrace=y^2(3-2y^2+3y^4) and W2=y^2(30-20y^2+3y^4)/30, with Atrace(0)=4 and A2(0)=1/30. Both weights are positive on0<y<1. Thus Atrace(q)>=4 and A2(q)>=1/30 for everyq>=0.

The trace reciprocal is cut-only. The shear reciprocal retains the exact S223 pole with R>0 and0.5798<r<0.5799. Substituting p=lambda^2+q into their complete first-sheet representations gives

1/Ftrace(lambda^2+q)=-integral rho_trace(tau)/(lambda^2+q+tau)dtau,
1/F2(lambda^2+q)=-R/(lambda^2+q+r m^2)
               -integral rho2(tau)/(lambda^2+q+tau)dtau.

The continuum densities in these expressions are the reciprocal densities, not the original growing stress cut. The shear pole frequency is now sqrt(q+r m^2); it is not deleted. These shifted scalar reciprocals have no right-half-plane Laplace pole. The scalar kernels need not have a uniform full half-lineL1 bound for the argument below.

Define the absolutely convergent primitives

Jtrace,q(t)=-integral rho_trace(tau)/(q+tau)(1-cos(sqrt(q+tau)t))dtau,
J2,q(t)=-R/(q+r m^2)(1-cos(sqrt(q+r m^2)t))
        -integral rho2(tau)/(q+tau)(1-cos(sqrt(q+tau)t))dtau.

The positive static measures have masses1/Atrace(q) and1/A2(q), respectively. Therefore

|Jtrace,q|<=1/2, |J2,q|<=60, Jtrace,q(0)=J2,q(0)=0

uniformly for allq>=0 and every fixed positivem. The integrals are continuous in time by dominated convergence. For Jq=diag(Jtrace,q,(3/8)J2,q),

||Jq(t)||<=45/2,
L[Jq](lambda)=diag(1/Ftrace(lambda^2+q),(3/8)/F2(lambda^2+q))/lambda.

Fubini is justified using the finite static measure, not by assuming absolute convergence of an oscillatory scalar kernel.

The forward factor is also well-defined as a causal distribution. Put c=1-y^2 and Omega_y,q=sqrt(q+4m^2/c). The kernels

G_i,q(t)=theta(t) integral W_i(y)sin(Omega_y,q t)/(c Omega_y,q)dy

are bounded and continuous, uniformly inq. Since Omega_y,q>=2m/sqrt(c), the elementary beta moments give

|G_2,q|<=9pi/(128m), |G_trace,q|<=27pi/(64m).

Their Laplace transforms are integral W_i/[4m^2+c(lambda^2+q)]dy. Thus the exact shifted forward distributions are
-A_i(0)delta-(D^2+q)G_i,q.
This keeps their otherwise divergent separated local/cut pieces together and automatically retains all initial distribution terms.

Independent60-digit tests use separately entered trace and shear cut formulas. Atq0,1,10000 and two real/complex Laplace frequencies, shifted radial inverses and their full spectral representations agree within1e-45. The recovered complete two-channel matrix product agrees with the identity within1e-38. These checks supplement, rather than replace, the analytic and positive-measure proof.
