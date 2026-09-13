# Full rolling-heavy block and its Euler contacts

Use the same logarithmic metric perturbation v, lapse n, shift divergence
b and Fourier conventions as S253/S254, now with hbar and mh nonzero.
The heavy normalized kinetic density before expansion is

    exp(3v) U(N+n)/(2(N+n))
        [mh+hdot-N^i partial_i h]^2,

with the entire potential (N+n)U(N+n)exp(3v)
[-n_mass(hbar+h)^2/2+JH(u,N+n)(hbar+h)/sqrt(kappa)] and the minimally
coupled spatial gradient. The raw full second variation is checked with
independent lapse jets through order2. The shift transport term is kept:
integral N^i partial_i h=-integral b h on a homogeneous background.

For both matter fields, c_i=Z m_i and w_i=Z_N m_i. The additional heavy
coefficients in the ENTIRE S254 four-mode block are

    d_H=partial_N[N U(-n_mass hbar+JH/sqrt(kappa))],
    vs_H=3N U(-n_mass hbar+JH/sqrt(kappa))
           -3(d_t+3H)c_H,
    mass_HH=-N U n_mass.

For M1, d_M1=0 and vs_M1=-3(d_t+3H)c_M1; the other matter masses vanish.
Replace the homogeneous base density in every lapse/volume coefficient
by the full L0 above. Then Cnn=L0_NN/2,
J=Cnn+3Theta^2/D-(w.w)/(2Z), Theta=-H D_N+B_N/2,
L0_lapse=3L0_N, and L2=4 partial_N(N C3) is unchanged.
Every resulting parameter is listed, including the whole time derivative
chain through N,H,mc,mh and hbar. This is a full simultaneous parameter
map into the complete S254 Hamiltonian, not a selected channel.

The raw-to-normalized heavy time boundary is
partial_t[kappa a^3 3c_H v h]. Together with the M1 and metric boundary,
it must be retained with the physical canonical map. The exact equations
imply the following lower contacts ONLY ON SHELL:

    L0_lapse=3E_N=0,
    Vvv=(3/2)E_volume=0,
    vs_i=3E_matter_i=0.

The metric identity uses cv=-18D H+9B and
E_volume=3L0-(d_t+3H)cv/3. All identities are derived before imposing
these equations. The lapse-heavy contact and rolling momenta are not
set to zero. The fast coefficient does not involve any of the four
vanishing lower contacts.

Holding hbar identically zero gives mh=mhdot=0 and recovers ALL S254
parameters exactly. In contrast, starting an unforced solution at
hbar=mh=0 generally gives mhdot=N^2 JH/sqrt(kappa), not zero. Its heavy
Euler equation cancels the volume contact. Conflating those two operations
would lose the source precisely where the on-shell argument needs it.
