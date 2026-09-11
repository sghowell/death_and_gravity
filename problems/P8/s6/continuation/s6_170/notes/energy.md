# Physical state-dependent energy and pressure differences

Let c_in and c_out be the SAME actual curved Hadamard states identified
in notes/hadamard.md. Each helicity covariance is a rank-one projector.
Their difference has trace norm2|beta_p|, unchanged by unitary evolution.

The physical T_hat00 bilinear is the original Hamiltonian with
omega(t)<=2E_global, E_global=sqrt(p^2+m0^2), because a>=1 and
M<=m+Delta<2m0. The physical volume factor is a^-3<=1.
Both helicities and N=42 copies therefore give

    |rho_in(t)-rho_out_state(t)|
       <=4N/pi^2 integral_0^infinity p^2 E_global |beta_p|dp.

Both states use exactly the same local subtraction and finite
counterterms; those terms cancel. No particular value of the
absolute curved vacuum stress or Newton-reference term is assumed.

For the compact transition allowance, E_global<=4Ec, p=4y and

    integral_0^infinity y^3/(y^2+m0^2)^(21/2)dy
         =2/(323m0^17).

Hence

    rho_compact<=262144N K^20(Delta+m tau)
                  /(323pi^2 tau^20 m0^17).

For the tail use Ap below p=m0 and B/p^5 above it. Below the split
E_global<=2m0; above it E_global<=2p. The complete radial integral
is at most

    A m0^5/2+2B/m0=4126mK^20/(15m0^17),

so

    rho_tail<=16504N mK^20/(15pi^2 m0^17).

Using pi^2>9 and the exact actual parameters gives a total below
1e-1090, with diagnostic value approximately2.6185804e-1095.
The ratio to the named kappa=1e800 scale is below1e-1890.
This is not a relative error against the vanishing bounce density.

The physical isotropic pressure operator before the mass rotation
is(q/3)sigma1. Its norm is q/3<=E_global/3, whereas the energy
allowance used2E_global. Therefore the pressure-state difference
is at most one sixth the same conservative energy bound.
Translation and rotation invariance give zero momentum density
and zero off-diagonal spatial stress.

For the external scalar mass, the appropriate difference identity is

    partial_t Delta rho+3H(Delta rho+Delta P)
       =Mdot Delta<bar psi psi>.

The fermion subsystem alone is not separately conserved. Neither
the source expectation nor absolute local curved stress is bounded
by this identity or by differentiating the energy inequality.
