# External theorem and a checked convention/state dictionary

The external input is [Gerard and Stoskopf's in/out-state
theorem](https://arxiv.org/pdf/2108.11955), version 2 dated
18 October 2021: Theorem 1.1, hypotheses H1-H3 in section
3.1.3, H4 in section 3.3, and the symbol definition in 4.2.
It provides pure Hadamard in/out states under its geometric,
mass-tail and asymptotic spectral-gap hypotheses. The theorem
is used, not re-proved or formalized by this repository.

Here the Cauchy surface is exactly Euclidean R^3. Its
injectivity radius is infinite and all curvature derivatives
vanish. The product spacetime has its trivial spin structure.
The lapse and its inverse are one, shift zero, and spatial
metric and inverse are the identity. All geometric asymptotic
differences are identically zero. Our mass-tail proof supplies
every required time derivative with symbol order -8.
Every spatial derivative of the mass is zero.

The two constant physical asymptotic masses are m +/- Delta,
strictly positive by the existing parameter enclosure.
Squaring either asymptotic Fourier Hamiltonian gives
(p^2+M_asym^2) times the identity, so zero is outside its
spectrum. The operator-norm projector limits just proved
identify our state with the one to which the theorem applies.

The paper uses signature (-,+,+,+). Our physical gamma
matrices gamma_+ have signature (+,-,-,-). Set, in the
source's orthonormal frame,

    Gamma_0=-i gamma_+^0,
    Gamma_j= i gamma_+^j,
    beta_source=gamma_+^0,
    m_source=-M(t).

Then the source Clifford relations hold,
Gamma_a^dagger beta_source=-beta_source Gamma_a,
and i beta_source Gamma_0=identity. More importantly,

    source slash + m_source
       = i gamma_+^mu partial_mu - M(t),

exactly our physical operator. With the paper's evolution
written as partial_t-i H_source, H_source=-H_physical.
Thus positive/negative spectral LABELS exchange. The
physical negative-energy occupied projector is not lost
by importing a convention-dependent sign. A negative
source mass here is a dictionary choice, not a negative
physical mass or a change of theory.

It follows that the specified free flat in/out states
are Hadamard. This conclusion uses the identified external
theorem and our hypotheses, not the previous finite-energy
bound. No numerical seminorm of the smooth remainder, local
stress prescription or interacting curved state follows
from this applicability check.
