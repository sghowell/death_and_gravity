# Actual full ADM Legendre reduction

Use the unchanged P8 +--- physical signature,

    ds^2=N^2dt^2-h_ij(dx^i+beta^i dt)(dx^j+beta^j dt),

with N>0 and positive h. Write v=sqrt(det h), F_ij=partial_i A_j-partial_j A_i and E_i=dot A_i-partial_i A0-beta^k F_ki. The ordinary connected Proca quadratic action is-F^2/4+m^2 A^2/2 in this signature. Direct contraction with the four-metric gives

    L=v h^ij E_i E_j/(2N)-Nv h^ik h^jl F_ij F_kl/4
       +m^2 v(A0-beta.A)^2/(2N)-Nm^2 v h^ij A_i A_j/2.

The conjugate momentum is Pi^i=v h^ij E_j/N. The Legendre transform contains Pi^i partial_i A0. Keep its boundary explicitly:

    Pi^i partial_i A0=-A0 divPi+div(A0 Pi).

After this spatial integration by parts, variation of the nondynamical A0 gives

    A0=beta.A-N divPi/(m^2 v).

Its substituted contribution is N(divPi)^2/(2m^2v)-beta.A divPi. The sign of the positive Schur term is checked directly. The full reduced Hamiltonian is

    H=integral N[h_ij Pi^iPi^j/(2v)
                 +v h^ik h^jl F_ij F_kl/4
                 +m^2 v h^ij A_iA_j/2
                 +(divPi)^2/(2m^2v)]
       +beta^i[Pi^jF_ij-A_i divPi].

There remain exactly the three physical Proca polarizations in the six-dimensional canonical phase space. Neither a fourth oscillator nor deletion of the longitudinal energy is allowed.

In three spatial dimensions write Bmag_i=epsilon_ijk F_jk/2. The cofactor identity for an arbitrary positive h gives

    v h^ik h^jl epsilon_ijr epsilon_kls/4
      =h_rs/(2v),

and hence magnetic density Bmag^t h Bmag/(2v). This follows from the two-epsilon determinant identity, not from assuming h diagonal. Literal arbitrary-SPD contractions independently test it.

Set N=1+n,h=a^2 exp(Q), Q any symmetric3x3 tensor, tau=trQ, B_Q=Q-tau I/2. Since v=a^3 exp(tau/2), all four lapse terms become

    H_N=N/2 integral[
      Pi^t exp(B_Q)Pi/a+curlA^t exp(B_Q)curlA/a
      +am^2 A^t exp(-B_Q)A
      +exp(-tau/2)(divPi)^2/(a^3m^2)].

For tracefree Q this exactly recovers S195. For Q=2psi I the factors are exp(-psi),exp(psi),exp(-3psi). In particular the constraint term is not independent of a spatial trace or lapse perturbation.
