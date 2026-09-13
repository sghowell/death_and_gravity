# Entire minimal-scalar stress pair and absolute cut normalization

Use eta=diag(+1,-1,-1,-1), canonical scalar mass squared n and two on-shell covectors k,l. With the universal1/sqrt(2E) on each mode reserved for invariant phase space, the symmetric Wick amplitude is

T_ab(k,l)=-(k_a l_b+l_a k_b)/2+eta_ab(k.l+n)/2.

There is no Wick exchange factor in this amplitude. The full matrix element has the expected exchange; in the connected two-current cut the single factor2 is included in the measure below.

At COM k=(E,0,0,kappa_r), l=(E,0,0,-kappa_r), E²=kappa_r²+n, the entire time row vanishes and

T_ij=k_i k_j-E² delta_ij.

For arbitrary homogeneous spatial Q, the exact scalar ADM Hamiltonian feature M_Q has momentum coefficient -trQ/2, gradient block (trQ)I/2-Q and mass coefficient trQ/2. Direct mode contraction gives
F_k^T M_Q F_l=-tr(Q T_spatial).
After normalization its pair is[ E² trQ-k.Q.k ]/(2E).
For Q=2vI this is v(2kappa_r²+3n)/E. Both the kinetic trace term and the mass term remain. A conformally improved scalar would give a different trace and is not this field.

Conservation of the full covariant pair follows by contraction with k+l, using k²=l²=n. Independent rational noncollinear Lorentz boosts test all components, not just a COM trace.

Rotational decomposition gives spin0 and spin2 eigenvalues

b=(tr T)²/3=(2kappa_r²+3n)²/3,
a=[tr(T²)-(tr T)²/3]/5=2kappa_r^4/15.

With s=4E² these are b=(s+2n)²/12 and a=(s-4n)²/120. Independent full9-by9 angular integration reconstructs the complete tensor a Pi2+b Pi0 for multiple massive ratios; it is not a tracefree interpolation.

The invariant two-body phase-space integral is beta/(8pi), beta=sqrt(1-4n/s). Connected Wick contraction contributes2. Converting the positive-energy Wightman cut to the Kallen-Lehmann density divides by2pi; the metric-current normalization contributes another1/4. Thus the final current cut is beta(a Pi2+b Pi0)/(32pi²). This yields exactly

rho_H2=beta^5 s²/(3840pi²),
rho_H0=beta(s+2n)²/(384pi²).

Both are zero through threshold. The positive retarded convention gives Im D_i(s+i0)=pi rho_Hi(s), and the three-subtracted function is D_i(z)=z³ int rho_Hi(s)/[s³(s-z)]ds. For unit-Frobenius shear its fourth-order reference factor is its independently fixed constant plus64pi² D_2(-p)/p². The Q=2vI trace normalization is12, hence768pi² for the scalar channel. The substitution s=4n/(1-y²) gives the full W_H2=y^6/30 and W_Htrace=y²(3-y²)². A cut fixes no finite subtraction coefficient by itself.

As a primary convention check, [Martin and Verdaguer, gr-qc/0001098v1](https://arxiv.org/abs/gr-qc/0001098v1), Eq5.8, gives the free scalar tensor noise spectrum. Setting minimal curvature coupling and converting their noise normalization to this current yields the same two densities. Their signature, subtraction convention and full gravitational propagator are not adopted. The normalization and all new inverse estimates above are derived here.

The cuts in this note are free external-metric Gaussian cuts, not interacting light/heavy scattering data or a UV completion.
