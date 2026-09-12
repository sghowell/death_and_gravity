# Original sharp mask, finite subtraction and all-transfer error

For the original pair mask chi_K(k)chi_K(P-k), K>=2m, S209-S211 give

    U_pair,K(P)-U_pair,K(0)
      =K^3 A3(P)+K A1(P)
        +F2(P)(K^2-m^2)/2+F4(P)log(K/m)+E_K(P),

with the complete actual curved A3,A1, the actual one-ball F2,F4, and

    |E_K[D,Gamma]| <1e40||D||L2 X46[Gamma]/K.

This error is uniform in all external P, not merely at fixed transfer. Complete conversion K2 and K0 artifacts cancel. That cancellation does not remove the physical F2 quadratic term or its finite lower-band term. The isolated leading-slot finite artifact is not re-added. Every finite odd endpoint and the evanescent odd UV coefficient remain accounted for.

Set

    S_K=K^3 A3+K A1+F2(K^2-m^2)/2+F4log(K/m)-Ffinite,
    Jren,K^anchor=H0+[Jactual,K(P)-Jactual,K(0)]-S_K.

The full one-leg contact cancels exactly in the bracket. This is an approximation anchored to the exact renormalized homogeneous current H0, not a claim about the rate for a separately unanchored homogeneous cutoff subtraction.

The exact error is

    Jren,K^anchor-Jren
      =[Known_K(P)-Known(P)]
        -[Known_K(0)-Known(0)]+E_K.

The same P0 specialization and source norm embeddings as in lift.md give each Known tail below1e54 M Z136/K. The third term is below1e40 M Z136/K. Therefore

    |Jren,K^anchor-Jren| <3e54 M[D]Z136[Gamma]/K.

The canonical display is12e-746/K. This establishes removal of the original auxiliary two-leg regulator in this anchored sector. It neither replaces that mask with a convenient physical regulator nor establishes a physical EFT cutoff, global-time bound, finite-gravity IR/Regge estimate or omitted-loop error.
