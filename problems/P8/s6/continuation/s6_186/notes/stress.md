# Actual local Proca stress and reference source

All fields and tests below use the physical orthonormal CD
frame. Let E_i=F_0i, B_i=(1/2)epsilon_ijk F_jk and
V_a=m A_a. Collect the ten real fields as Z=(E,B,V).

For signature(+---) the actual canonical action gives

    T_ab=-F_ac F_b^c+eta_ab F_cd F^cd/4
         +V_a V_b-eta_ab V_c V^c/2.

The complete components are
T00=(E^2+B^2+V0^2+Vsp^2)/2,
T0i=(E cross B)_i+V0 Vi, and
Tij=-Ei Ej-Bi Bj+Vi Vj
    +delta_ij(E^2+B^2+V0^2-Vsp^2)/2.
The exact covariant tensor and these component formulas
are independently reconstructed. Write Tab=Z^T M_ab Z.
Every real symmetric ten-by-ten M_ab has row/operator
norm at most1/2. Thus for any complex tensor f,

    ||M_f||op<=sum_ab|f_ab|/2<=2||f||F,

including all16 entries, not just the diagonal currents.

For k along the third axis, one transverse mode has
phase-space readouts (E1=p, B2=i(k/a)f, V1=mf).
The longitudinal mode has
E3=(m/omega)p, V0=-i(k/a)/omega p, V3=omega f.
Every readout has a common a^-3/2 prefactor. The second
transverse mode is its orthogonal rotation. These follow
from gT^2=a, gL^2=a m^2/omega^2 and the actual temporal
constraint A0=-div pi/(a^3m^2), with p=f'-d f.

The complete one-particle energies are
(|p|^2+omega^2|f|^2)/(2a^3) in each polarization.
No longitudinal magnetic mode or independent A0 oscillator
is introduced. Arbitrary momentum directions only rotate
the physical components, preserving the norm bounds.

The local stress smearing uses a^3 dt dx. This volume
cancels both a^-3/2 readout prefactors exactly in each
pair amplitude, before differentiating the test in time.

At the reference clock the FULL source has both R-1=0
and Q=0. Consequently its first metric or scalar
variation is zero. The source contact and linear
centered-vector coupling supply no extra first-force
noise there. This leaves the ordinary quadratic
Proca stress, not zero metric noise.
