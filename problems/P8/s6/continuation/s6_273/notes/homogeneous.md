# Absolute canonical variables and every homogeneous force

Let alpha=log(a), p=PV/(3a^3), pm=PM/a^3, ph=PH/a^3 and eta=10^100 h.
The normalized homogeneous Hamiltonian density f is the restriction of
the FULL Hbar to homogeneous isotropic fields, after the complete
temporal elimination. Keep R,F,j,I,Iu, all mass terms and Hclock.
Use the actual root f_N=0; the full lapse envelope theorem is essential.

The original canonical one-form gives

    alpha'=f_p/3,
    p'=-f+pm f_pm+ph f_ph,
    pm'=-pm f_p,
    eta'=10^100 f_ph,
    ph'=-10^100 f_eta-ph f_p,
    M1'=f_pm.

All six equations are independently derived from a^3 f in the canonical
variables before imposing the lapse root. In particular the dilution
terms and the physical heavy factor10^100 are not omitted.

At the actual prescribed reference p=-2Hclock, pm=ell, eta=ph=0, N=1,
the clock rates agree with the reference derivatives except

    p'_clock-p'_reference = -PRESSURE_fixed(u).

The full clock constraint also has the retained residual
3 PRESSURE_fixed/[2(1+u^2)^3]-RHO_fixed. The reference is therefore not
declared to be a solved classical background. Its true lapse displacement
is bounded using the complete constraint. Fixed profiles are not varied
as evolving quantum means.

## Canonical split, not a missing moving-chart term

The full S267 one-form is evaluated with absolute homogeneous variables
and zero-mean INPUT fluctuations:

    gamma=a^2 exp(2v)Q^-1,  Pi_v_total=PV+deltaPi_v,
    M1_total=M1_hom+deltaM1, Pi_M_total=PM+deltaPi_M,
    h_total=eta/10^100+delta h, Pi_H_total=PH+deltaPi_H.

For example, integral(PV+deltaPi_v)d(alpha+v)
equals Vol PV d(alpha)+integral deltaPi_v dv. Its cross terms vanish
because the input fluctuations and their variations have zero mean.
The complete nonlinear cotangent lift already preserves this one-form;
generated means and harmonics in reconstructed fields are not discarded.

The same statement holds for the matter channels and full shape/vector
pairs. Therefore the hybrid canonical form has no extra Y'-dependent
quantum connection. At the old external reference, the difference between
the absolute Hamiltonian and the old moving-chart Hamiltonian is the
known scalar connection kappa Vol(Hclock PV+ell PM). The zero-mean
fluctuation contacts vanish only AFTER restriction. The original primitive
and any exact boundary are counted once.

## Quantum-force pullback

Let k_i=<partial_i K_Y> for the complete CENTERED bounded quantum
generator, and Dhom=kappa Vol a^3. The exact additional rates are

    delta alpha' = k_p/(3Dhom),
    delta p' = (-k_alpha/3+pm k_pm+ph k_ph)/Dhom,
    delta pm' = -pm k_p/Dhom,
    delta eta' = 10^100 k_ph/Dhom,
    delta ph' = (-10^100 k_eta-ph k_p)/Dhom,
    delta M1' = k_pm/Dhom.

They follow by changing variables from the canonical derivatives of the
full expectation energy. The generic functional check uses an independent
alpha argument until AFTER partial differentiation; it does not conflate
a partial alpha derivative with the total derivative of density momenta.
All six identities are tested. The full M1 charge a^3 pm is exactly
conserved, including these quantum corrections.

The scalar center h0 remains the complete classical homogeneous energy.
Removing its phase for a comparison norm does not remove its force.
A parameter-dependent scalar has a generally nonzero classical derivative.

## Uniform rate estimates

For the full classical rates V_i, differentiate first at fixed N.
At the implicit root,

    partial_z V_i = (V_i)_z-(V_i)_N C_z/C_N.

The complete source bounds, |CN|>2 and every retained derivative give
the full N-force and five-dimensional row-sum Lipschitz ceilings1e112.
The alpha column is zero for the normalized homogeneous density equations;
M1 is cyclic and is reconstructed afterward. The large heavy frequency
is retained in these ceilings, not replaced by a mass-independent claim.
The source-controlled reference lapse shift then gives a full reference
rate error below1e-280.
