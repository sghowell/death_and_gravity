# Exact root-free local action

Keep the complete constrained Hamiltonian M=diag(V,K), with K positive
symmetric and V=omega^2 K^-1. The scalar omega is positive. These two
independent variables cover this instantaneous Hamiltonian class; no
commuting polarization ansatz is introduced.

For B=K^(1/2), L=B^-1 B', R=(L^t-L)/2,
P=(L+L^t)/2 and S=rho I-P, rho=omega'/(2omega),
differentiate B^2=K in the actual order. One obtains

    P=B^-1 K' B^-1/2,
    D P=B^-1 K'' B^-1/2-2P^2, D=partial_t-ad_R.

Let A=K'K^-1 and p=omega'/omega. Similarity, not commutation, gives

    B S B^-1=s=pI/2-A/2,
    B(D S-pS)B^-1=t=s'-p s
       =(p'-p^2)I/2-A'/2+pA/2.

Thus the complete local marker terms can be written without a matrix root:

    L0=-d omega/2,
    L2=tr(s^2)/(4omega),
    L4=tr(t^2+s^4)/(16omega^3).

Their interpretation as the local action for the actual fourth-order
current is proved by the independent variational calculation in
variation.md. It is not inferred from a vacuum-phase heuristic.

The ordered Riccati recurrence also gives
r1=iS/(2omega), r2=(DS-pS)/(4omega^2) and

    r3=i[S^3-D^2S+3pDS+(p'-2p^2)S]/(8omega^3).

The trace phase expression at these orders is consistent with L2,L4
after compact integration by parts, but that consistency is not used
in place of the full independent current variations. A constant squeeze
checks the sign against -sqrt(omega^2-S^2)/2 in one mode.

Only finite local orders0,2,4 are represented here. This is not a
convergent infinite adiabatic series, a complete determinant, a new
state, or a scalar replacement for the constrained vector.
