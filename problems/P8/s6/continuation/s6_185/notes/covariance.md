# Same all-order CD state: explicit one-particle covariance

Keep the fixed physical CD metric a=(1+t^2)^2 on I=[-1/2,1/2],
the canonical mass1000 and the actual S6.55 all-order Cauchy
state already matched in S6.176. A scalar-history change at
fixed metric changes the coherent source mean but not this
connected covariance. No flat vacuum or finite-WKB state is
substituted.

For each transverse or longitudinal canonical mode, the
source-pinned W8 comparison satisfies
omega/2<W8<3omega/2 and |W8'/W8|<4.
The actual physical canonical rate is
dT=H/2, dL=(1/2+z)H, z=q/(q+m^2); hence |d|<3.
The reference mode and physical momentum are
fref=exp(-i integral W8)/sqrt(2W8),
pref=(-iW8-W8'/(2W8)-d)fref.
Consequently |fref|^2<=1/omega and
|pref|^2<3omega, using m>=1000.

The unchanged all-order state's exact mixing coefficients
obey |beta|<1 and |alpha|^2-|beta|^2=1. Their transport
constraint eliminates derivatives of alpha,beta in the
mode readout, so both f and p are the same linear
combinations of reference plus/minus modes. Therefore
|alpha|+|beta|<3 gives

    |f|^2<9/omega, |p|^2<27omega

at every time and every spatial momentum. This uses the
full state and its evolved preparation envelope, not a
vacuum reset or a truncated high-momentum expansion.

Let k denote comoving momentum, q=|k|^2/a^2 and
omega=sqrt(m^2+q). The actual canonical maps are
gT^2=a, gL^2=a m^2/omega^2, A_sp=f/g and pi_sp=g p.
The free temporal constraint is A0=-div pi/(a^3m^2).
In an orthonormal physical frame, the trace of the
complete three-mode covariance matrix at each time
is bounded by

    a^-3 [27/omega+36q/(m^2 omega)].

The two transverse terms give18/(a^3 omega); the
longitudinal spatial term gives9omega/(a^3m^2), and
the temporal term gives27q/(a^3m^2omega).
The nonzero temporal and longitudinal components are
essential. Cross components are controlled by positivity
and the Euclidean trace bound.

For a smooth compact test four-vector f in physical
orthonormal components, the field smearing uses volume
a^3 dt dx. Cauchy-Schwarz in time, the exact positive
mode covariance, and spatial Plancherel yield

    C_A(f,f)<=T Amax^3 [
        27/m ||f||L2(dt dx)^2
        +36/m^3 ||grad_spatial f||L2(dt dx)^2 ],

where T<=1 and Amax=25/16. Since Amax^3<4, one may use
coefficients108/m and144/m^3. This is a full all-momentum
state-dependent UPPER covariance bound. It is not a QEI
lower bound, a flat spectral formula or stress-tensor
noise. Metric index flips preserve the Euclidean norm;
proper-frame conversion of lower spatial test components
only decreases their coordinate norm since a>=1.
