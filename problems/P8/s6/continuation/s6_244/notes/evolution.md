# Full scalar evolution and current contacts

Let v=a³exp(tr Q/2), K=1/v, V=v omega² and omega²=n+a^-2 k^t exp(-Q)k. With canonical x=(field,pi), J=[[0,1],[-1,0]], the complete Hamiltonian matrix is M=diag(omega²/K,K). There is one scalar oscillator per momentum, not a vector polarization sum.

The symplectic time-dependent frame T=diag(sqrt(K/omega),sqrt(omega/K)) gives

T^-1 J M T-T^-1 T'=[[s,omega],[-omega,-s]],
s=(lambda+3H+tr Q'/2)/2, lambda=omega'/omega.

The physical annihilation solution uses -iomega. Its normalized pure graph r=x+iy satisfies r'=2iomega r+s(1-r²). For d=1-x²-y² the ENTIRE balanced real covariance is

Sigma=[[1+2x+x²+y²,-2y],[-2y,1-2x+x²+y²]]/(2d).

The code checks all four entries of Sigma'=A Sigma+Sigma A^t and det Sigma=1/4. In particular Sigma_QP=-y/d is retained. Exact smooth real Hamiltonian evolution preserves positivity/purity, hence |r|<1; the quantitative proof below improves this to |r|<1/10. No finite adiabatic polynomial is evolved as the physical covariance.

For a logarithmic spatial variation D, delta K/K=-tr D/2 and delta omega² is the full Frechet derivative of exp(-Q). Thus

N_D=T^t M_D T/omega=diag(tr D/2+delta_D omega²/omega²,-tr D/2),
J_D=-omega tr(N_D Sigma)/2=-tr(M_D C)/2, C=T Sigma T^t.

No tracefree restriction is imposed. The matrix exponential derivative retains noncommuting D and Q.

With the initial Cauchy covariance fixed, the exact canonical tangent is

C1'=A C1+C1 A^t+A1 C+C A1^t, C1(t0)=0,
C2'=A C2+C2 A^t+2(A1 C1+C1 A1^t)+A2 C+C A2^t, C2(t0)=0.

Therefore the complete first current derivative is -tr(M_D C1+M_DG C)/2. The second is -tr(M_D C2+2M_DG C1+M_DGG C)/2. All instantaneous contacts and factors of two are present. A separately evolved covariance IVP tests both derivatives against the full noncommuting matrix-exponential Hamiltonian, with contacts nonzero at readout.

For the finite reference set rhat=sum from j1 to10 of r_j, with

r1=-s/(2iomega),
r_(n+1)=[r_n'+s sum_(j+l=n)r_j r_l]/(2iomega).

The complete residual is F=r10'+s sum_(j+l>=10)r_j r_l, where both indices range1..10. With e=r-rhat, direct differentiation gives the exact first and second error equations in the report. The second equation retains -2s e1². It is not a linearization of the error equation.

For real histories r_j has phase i^j times a real coefficient. Marker reversal interchanges rhat and its formal conjugate. The DIAGONAL spatial current is even in the marker. The off-diagonal covariance is odd, generally nonzero, and is never declared even.
