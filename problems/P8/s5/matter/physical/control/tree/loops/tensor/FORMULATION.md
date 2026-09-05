# P8 S5.10.CD — Constant-Weyl tensor matching candidate

Adopted 2026-09-05. The separately named candidate is
S_CD/M1+cC integral sqrt(-g) C², with constant cC and the original physical
matter metric. It does not overwrite the old covariant action or any prior
certificate. It pins/replays S5.9 SHA256
`da7e0b445245dd87369380c087dc8a6a45361b68e2b4711f2c6ffa989db42fb1`.

## Scope and physical variable

The full first variation of integral C² is zero on flat FLRW, so the old
CD/M1 background remains an exact solution of this new action. Only its
linear tensor sector is investigated here, separately for both unit-norm
TT polarizations. Scalar/matter constraints and new nonlinear interactions
are not certified.

P8 has signature +---, Einstein term -M²R/2 and unchanged canonical-free
matter. For the physical metric perturbation gamma and conformal time,

    S_T = (M²/8) integral a²(gamma'²-k² gamma²)
          + (cC/2) integral (gamma''+k² gamma)².

Let beta=cC/M², q=k²/a² and B=Hddot-2H Hdot. The selected first-order
order-reduced representative is

    gamma_ddot+(3H+8 beta B)gamma_dot+(1-16 beta Hdot)q gamma=0.

Its exact two-data solutions match the formal perturbative tensor equation
through first order in beta and have an explicitly bounded second-order
residual in the original fourth-order equation. This is not an existence,
uniqueness or approximation theorem for exact four-data higher-derivative
branches. Solving this particular reduced equation exactly is an explicitly
recorded finite-order representative, not a preferred all-orders completion.

Action reduction uses the off-shell map
gamma=y+2 beta E0[y]-8 beta H y_dot+O(beta²). Only on the perturbative
branch E0[y]=O(beta) can its E0 term be dropped. The reduced y action has
kinetic coefficient 1+16 beta H²; its speed coefficient is not the one
for the original physical gamma. The observable map must be retained.

## Centre-data finite-window contract

For every finite centre t0 let ell0=tau sqrt(1+u0²), a(t0)=1 and
|t-t0|<=ell0/100. Require centre momenta
10^11<=k_com ell0<=4*10^12. This covers the bounded tensor momenta used
by S5.8; it is neither a loop regulator nor a global fixed-comoving band.
Supply independently the constant-coefficient bound

    epsilon = |cC|/(M tau)² <= 10^-30.

No finite value of cC is inferred from the S5.9 pole or logarithmic running.
The compared solutions have the same physical gamma(t0), gamma_dot(t0).
Equivalently they have the same old canonical Y,P at the centre, where

    Y=(M/2)a^(3/2)gamma, sigma=(t-t0)/ell0, P=dY/dsigma,
    w0=ell0²[q-(3/2)Hdot-(9/4)H²],
    E0=(|P|²+w0|Y|²)/2, N_t=sqrt(2 E0).

No equality of redefined y data, corrected canonical momenta or quantum
vacuum states is assumed. The rigorous centre-relative results are

    3/5 <= E0_reduced(t)/E0_reduced(t0) <= 5/3,
    N_t(z_reduced-z_free) <= 10^15 epsilon N_t0(z_initial),
    N_t(z_reduced-z_free-epsilon0 z_first)
        <= 10^30 epsilon² N_t0(z_initial),

where epsilon0=cC/(M ell0)² is signed and z_first is the derivative of the
chosen representative's solution with respect to epsilon0 at zero, with
zero correction to the centre data. The energy ratio is for nonzero data;
the norm inequalities include zero data. Frequency enhancement is retained.

The original-equation residual satisfies

    ell0² |(M/2)a^(3/2) Efull[gamma_reduced]|
        <= 10^38 epsilon² N_t(z_reduced).

It is a residual bound, not a bound on a singular fourth-order inverse.
At M tau=10^324 and the separate assumption |cC|<=1, the three bounds are
respectively 10^-633, 10^-1266 and 10^-1258. These are sufficient, very
loose conditional estimates, not a naturalness prior or optimal cutoff.

## Causality and quantum limitations

Hdot changes sign at |u|=1. Thus every nonzero constant cC gives both signs
of the leading physical two-derivative speed shift somewhere on CD. The
coefficient is positive under the chosen smallness condition, but one sign
is above the old unit coefficient. This is not a full quantum or UV
superluminality verdict. A reduced low-energy coefficient is not the
high-frequency front velocity, and the additional roots of the resummed
fourth-order truncation are not established physical EFT states.

The candidate is only one local operator, not the isolated matter loop's
complete effective action. Finite matching, other curvature operators,
nonlocal/anomaly/state terms, scalar and nonlinear effects, unknown higher
orders, backreaction and UV completion remain open. Constant-scale running
cannot be substituted for a finite matching coefficient, and mu(t) is not
inserted into this constant-cC action. Old certificates retain their old
actions and scopes.
