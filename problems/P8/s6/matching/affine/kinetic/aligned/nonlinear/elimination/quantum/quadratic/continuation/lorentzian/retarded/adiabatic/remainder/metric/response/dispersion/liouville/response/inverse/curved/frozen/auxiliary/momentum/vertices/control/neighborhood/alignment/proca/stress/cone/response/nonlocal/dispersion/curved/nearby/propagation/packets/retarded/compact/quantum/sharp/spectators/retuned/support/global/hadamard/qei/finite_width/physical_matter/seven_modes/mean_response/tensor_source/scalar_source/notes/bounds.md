# Whole-strip coefficient, state and mean bounds

Take |u|<=1/100 and 1<=k<=2, with the smooth state
addition supported strictly inside the corresponding band.
Every generator/source/observable entry is rational in u,k.
Its denominator has a positive constant and only positive
even-time monomials after exact cancellation. The numerator
absolute-coefficient sum at |u|=1/100, |k|=2 therefore gives
a rigorous upper bound divided by that denominator constant.
The checker rejects expressions without this proof structure.
No momentum/time grid or floating fit is used.

The actual natural scalar generator row sums are below
1,2,87,5 respectively, so its infinity norm is below 90.
Its fundamental matrix on either half of the strip has
norm at most exp(90/100)<exp(1)<3; the factorial series
gives exp(1)<3. The added covariance has every entry bounded
by 9eta after band integration: each row Euclidean norm
is bounded by its row sum, and the anchor covariance is I4.

The half entry sums of the four source Hessians are below
38,1,41,6. Hence their expectations are bounded by

    |F_N|<350eta, |F_xi|<10eta,
    |F_p|<400eta, |F_psi|<60eta.

The unchanged classical mean matrix has infinity norm
below one on this strip. Using J>3/2, |alpha|<1/10,
ell<=1/10, |beta|<=1/20 gives complete forcing below
410eta and, with zero anchor scale/trace,

    max(|xi|,|dp|)<=410eta/99<5eta.

The constraint and matter equation give

    |n|<120eta, |xi+n/(2h)|<66eta,
    |delta psi'|<70eta, |delta psi|<eta.

The intrinsic density and pressure Hessian half entry sums
are below 4 and 3. Their state expectations have magnitude
below 36eta; the additional mean contribution is below
(3/100)66eta<2eta. Thus both complete leading physical
matter observable differences have magnitude below 40eta
in units M^2/tau^2. These bounds concern means and quadratic
expectations, not every configuration of a Gaussian field.

All rational entry estimates and every rounded inequality
are independently recomputed. The displayed constants are
nonoptimal and apply only to this fixed compact band and
time strip. The regular all-finite-time continuation of the
linear mean system does not extend these numerical bounds
to the infinite tails.
