# Physical lapse observable, canonical normalization and local pole

Retain the S6.42 physical metric, light-field dictionary and vector
action. Work first in units tau=1, with the canonically normalized
vector Wc=sqrt(zeta_physical) W and m=m0*tau. Write a_s for the
physical scale factor, and a_m(N),b_m(N) for the two mass coefficients.
They equal one on the clock but their lapse derivatives are

    alpha=4/(9h), beta=28/(81h), h=(1+u^2)^3.

This is the actual S6.47/48 dictionary, not freely chosen Proca
couplings. The effective vector source and its first variation vanish
on this clock. Its Gaussian mean is zero. The field translation
W=T-B(phi,X)dphi has unit translation Jacobian; changing N at fixed
T adds a term linear in W to the first variation, whose Gaussian
expectation vanishes. The following quadratic expectation is thus
the retained vector's contribution to the physical lapse equation.
This statement does not discard auxiliary contact/measure terms of
the full affine path integral or supply other field determinants.

Let q=k_com^2/a_s^2, omega^2=q+m^2, Wi=i*k_i*sigma in the
longitudinal sector. At fixed physical spatial metric the per-mode
real quadratic Lagrangians are

    LT=a_s*w'^2/(2N)-N*a_s*(q+m^2*b_m)*w^2/2,
    LL=a_s^3*[q*(sigma'-W0)^2/(2N)
              +m^2*a_m*W0^2/(2N)-N*m^2*b_m*q*sigma^2/2].

The temporal equation gives W0=q*sigma'/(q+m^2*a_m). The energy
is -a_s^-3 times the N derivative before N=1. Taking that derivative
before eliminating W0 or after eliminating it gives exactly the same
answer, as verified without using a dynamical W0 equation. Only the
first mass jets are needed in this first variation. In particular
setting a_m=b_m=1 before variation is incorrect for this observable.

The physical canonical maps and rates are

    w=v/sqrt(a_s), sigma=v/sqrt(a_s^3*m^2*q/omega^2),
    z=q/omega^2, dT=H/2, dL=H*(1/2+z).

Polarization vectors are unit normalized. Complex normalized modes
replace real squares by absolute squares. With two transverse and
one longitudinal polarization, the densities per polarization are

    rhoT=[|v'-dT*v|^2+(omega^2+beta*m^2)|v|^2]/(2*a_s^3),
    rhoL=[(1-alpha*z)|v'-dL*v|^2+omega^2*(1+beta)|v|^2]/(2*a_s^3),
    pT=[|v'-dT*v|^2+(q-m^2)|v|^2]/(6*a_s^3),
    pL=[(1+2z)|v'-dL*v|^2-omega^2*|v|^2]/(6*a_s^3).

Pressure is the isotropic spatial variation (3Na_s^2)^-1 dL/da_s
at fixed k_com, not fixed q. It has no alpha,beta term since X is
unchanged in this spatial variation. The zero-momentum canonical
limits of both rho and p agree between sectors. The sigma chart is
not used at k_com=0; the three homogeneous coordinate-vector modes
provide its regular replacement. In the alpha=beta=0 control the
ordinary Proca expressions are recovered. A useful primary comparison
is [Maranon-Gonzalez and Navarro-Salas, Eqs. (46)-(49)](https://arxiv.org/html/2412.01963v1),
with their mode Wronskian 2i converted to ours i. Their ordinary Proca
energy does not include the present clock mass variation.

On I=[-1/2,1/2], 0<=z<=1 and h>=1. Thus the two positive energy
weights A,B, defined by rho=(A|v'-dv|^2+B|v|^2)/(2a_s^3), satisfy

    5/9<=A<=1, omega^2<=B<=109*omega^2/81<3*omega^2/2.

Replacing the signed pressure weights by their absolute values gives
at most 9/5 times this energy form: (1+2z)/3<=1<=(9/5)A
and omega^2/3<=B/3 in the longitudinal sector; transverse is smaller.
The vector component need not separately obey the ordinary conserved
Proca stress identity: clock-dependent masses exchange with the light
equation. No stress-conservation or full gravitational equation is
inferred from these forms alone.

## Independent zero-derivative pole check

In spatial dimension 3-2epsilon_DR, the scalar zero-point integral is

    I = mu^(2epsilon_DR)*(m^2)^(2-epsilon_DR)
        *Gamma(epsilon_DR-2)
        /[2*(4*pi)^(3/2-epsilon_DR)*Gamma(-1/2)].

The direct Gamma residue is -m^4/(64*pi^2*epsilon_DR); the harmless
mu factor does not affect it. At derivative order zero the two energy
forms integrate, per polarization, to

    rhoT=I+beta*m^2*I_mass_derivative,
    rhoL=I+(beta-alpha)*I/2+alpha*m^2*I_mass_derivative.

Because the pole of I_mass_derivative is 2I_pole/m^2, their sum has
relative pole 3+3alpha/2+9beta/2=3+20/(9h). This exactly equals
C+C_N from the frozen local determinant potential: the lapse energy
of a local term V(N) is V+NV_N. Continuing the polarization count to
d dimensions can change its finite coefficient but not this pole.
Omitting the clock jets misses 20/(9h). This comparison neither
derives a finite subtraction scheme nor identifies all derivative
counterterms in the physical energy equation.
