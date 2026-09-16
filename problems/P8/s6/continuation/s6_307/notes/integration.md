# Uniform total variation of the one-graviton difference

DefineTV(x) by integrating
|rho sum_pol|M5/A_B|^2-sum_pol|S0|^2|
over all emitted directions and0<omega<=x. The angle-integrated
one-graviton measure isomega*domega/(4pi^2). Taking the pointwise
absolute value bounds the absolute signed finite rate.

Split ata=sqrt(delta)/192 whena<x. The two regions cover every
nonforward angle without a fixed lower bound ondelta.

In the low band putB=2e8. The full remainder is at mostB/sqrt(kappa),
the common soft current at most64/(sqrt(kappa)omega),rho<=1 and
|rho-1|<=2omega. After both polarizations the integrand is bounded by

[(256B+16384)/omega+2B^2]/kappa.

Multiplying the measure and integrating up to min(x,a), bounded byx,
gives[(256B+16384)x+B^2*x^2]/(4pi^2*kappa).
Sincex<=1/8 and4pi^2>36 its coefficient ofx is below
1250012800004096/9<2e14.

For the high band use the direct gravity amplitude bound, not the
low-energy remainder. S304 gives|S0|<=160sqrt(delta)/(sqrt(kappa)omega).
The positive full-Born weights and the matter remainderBm=330000 imply
the pointwise two-polarization absolute difference is at most

[6*160^2*delta/omega^2+4Bm^2
 +2D^2*delta^2/(omega^2*(delta+omega^2))]/kappa,

whereD=1e17. This includes every matter-gravity interference by the
convex squared-amplitude inequality.

Use delta^2/(delta+omega^2)<=delta. Integratingomega*domega yields
an upper numerator

(6*160^2+2D^2)*delta*ln(x/a)+2Bm^2*x^2.

When this band is nonempty let y=sqrt(delta)<192x.
The elementary function y^2*ln(192x/y) is zero at both endpoint limits
and maximized aty=192x/sqrt(e), with value(192x)^2/(2e).
It is therefore strictly below18432x^2. Division by4pi^2>36 gives
the high coefficient10240000000000000000000000006128643200<2e37.

The result is

TV(x)<[2e14*x+2e37*x^2]/kappa<3e36*x/kappa.

In addition S304 independently givesTV(x)<1e32/kappa, so the minimum
of those two valid upper bounds is also valid. The new bound vanishes
uniformly asx->0 even whiledelta->0 by an arbitrary path.
It does not make the forward Born cross section finite or exchange
the already ordered virtual-real dimensional-regulator limit withx.
