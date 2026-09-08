# P8 S6.57 — Fixed-vector-stress homogeneous response

Keep the new S6.56 action, its original physical metric and the
S6.55 Hadamard vector state and finite prescription. This step
does not change the action or the quantum state. Original P8
remains open.

On I=[-1/2,1/2], treat the already defined actual vector energy
and pressure on the classical history as a fixed first-order
source. Normalize them by M²/tau² and their first u derivatives
by the corresponding physical M²/tau³ scale. Let eta0 bound
both values and eta1 both derivatives, with 0<=eta0,eta1<=10^-6.

The literal homogeneous action is restricted before eliminating
a nonzero spatial mode. Its physical stress forcing is

    F_v=3p, F_n=3delta*p-rho, delta=1/(2h).

The second term includes the physical spatial-metric response to
the lapse. The compatible clock source is the unchanged Ward
source from S6.54--S6.55, not a separately adjustable datum.

With zero independent canonical perturbations at u=-1/2,
there is a unique regular first-order homogeneous response.
The lapse is solved, not also imposed to vanish initially.
A two-component normalized phase system has operator norm
below 3/2 and forcing below 49eta0 over the whole interval.

Continuous bounds give

    |v|, |4a³(p_v+3ell*s)| <=131eta0,
    |n| <=94eta0, |n'| <=48000eta0+4eta1,
    |s| <=44eta0.

The exact positive physical point-chart lift of this first-order
response has fractional scale-factor change <=450eta0 and
Hubble change <=(98000eta0+8eta1)/tau. Its lapse remains in
the original clock tube. Its Hubble rate retains opposite signs
at u=+-1/4, so this approximate metric has an interior local
minimum of the scale factor.

At the actual vector example L=M*tau=10^12,R=m0*tau=1000,
eta0,eta1<10^-14. Hence the lifted scale changes by less than
4.5*10^-12 and its Hubble rate by less than 10^-9/tau.

This solves the stated fixed-source **linear** response, not
the self-consistent quantum equations. The lift is a specified
higher-order completion of linear fields, not a nonlinear
solution. No quantum-feedback norm, complete residual bound,
nondegenerate/unique bounce, global quantum history, corrected
perturbation cone, cutoff or V/G/B gate is established.
