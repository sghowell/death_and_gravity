# Actual first formal stationary quantum response

Use S6.182's new fixed action, not an alteration of its
unchanged predecessor. Let r=rho_ref/kappa, p=P_ref/kappa
and S=r+p. The actual unretuned physical vector currents
are (-delta r-3r v, 3delta p+3p(n+3v)).
The Hessian of the new FIXED coefficient contributes

    [ -S, 3r; 3r, -9p ] (n,v).

Their sum is therefore

    Qstar(n,v)=(-delta r-S n, 3delta p+3S n).

The code checks this literal physical sum and independently
maps it to the complete prepared adapted response in
S6.182. It retains the same nonlocal scale channel and
all contacts. Reference r,p are fixed functions, not
re-evaluated under variation to cancel the new response.

The source-pinned physical conditional stress bound from
S6.179 is less than5e-795 for each component per unit
joint physical C10 norm. The reference S6.176 stress
bounds give |S|<2e-770. Consequently

    ||Qstar h||C0 <7e-770 ||h||C10,
    ||Bcl Qstar h||C0 <4e-765 ||h||C10.

The complete homogeneous source-centered affine vector
mean force is identically zero; no omitted additional
classical mean term is silently placed in this Qstar.

Introduce a FORMAL alpha multiplying BOTH the fixed new
coefficient and the specified conditional Gaussian
functional. The retained reference mean is stationary
coefficientwise for every alpha, and T0 is the unchanged
classical reference. Formal coefficient matching in

    (T0+alpha Qstar)(h0+alpha h1+...)=g

gives T0 h0=g and T0 h1=-Qstar h0.
Both equations have unique smooth prepared solutions
by the actual classical propagator. Source-pinned
Hadamard mode response is smooth on smooth prepared
inputs; its numerical bound is only C10-to-C0.

The proven classical C10 estimate and the actual
composition bound give the sharper exact majorant

    ||h1||C0 <53000*(7e-770)*(1e96)||g||C10
             =3.71e-669||g||C10
             <4e-669||g||C10.

This is a well-defined formal coefficient, not a claim
that an exact solution branch is differentiable at
alpha=0. No existence of that branch is needed to define
h0 or h1. In particular this number does not bound the
finite-alpha=1 correction or the omitted second-order
term. Smallness C10-to-C0 is not a small endomorphism
on C0; no Neumann contraction or causal-pole deletion
is justified by this result.
