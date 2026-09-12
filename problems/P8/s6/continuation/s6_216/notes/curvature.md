# Full-volume curvature Hessians without a tracefree shortcut

Use conformal time eta, g=a^2[-deta^2+exp(Q)_ij dx^i dx^j], and retain exp(tr(Q)/2) in the volume. Detector and source have opposite spatial Fourier momenta. Move detector time derivatives onto the complete coefficient and source before converting d/deta=a d/dt and dividing the coordinate density by a.

Let AF=a^(d-2), PF=a^(d-4), H=a'/a and E=-T+2V+S-UD-UG. The spatial Einstein operator is AF p^2 E Gamma/2. Its background scalar is R0=d[2H'+(d+1)H^2]. The linear scalar curvature of a direction Q Gamma is

    tr(Q)[Gamma''+(d+1)H Gamma']
       +(p^2/a^2)[tr(Q)-eQe]Gamma.

Thus the full R^2 spatial Hessian is2R0 H_R plus the spatial difference of2a^d delta R_D delta R_G. Define AD=S-UD and AG=S-UG. Integration by parts gives

    H_R2=2R0 H_R
      +2AF p^2[(AD+AG)Gamma''
        +((d-5)AD+(d+1)AG)H Gamma'
        -3AD(H'+(d-2)H^2)Gamma]
      +2PF p^4(S-UD-UG+W)Gamma.

For Ricci^2 a direct conformal ADM decomposition avoids a tracefree assumption. With L=a_eta/a and K=exp(-Q)(exp Q)_eta/2,

    (Ric_g)^i_j=R3^i_j+K'^i_j+tr(K)K^i_j+(d-1)L K^i_j
                  +[L'+(d-1)L^2+L tr(K)]delta^i_j.

Ric00 has no spatial derivative. Ric0i at first order is(div Q'_i-grad_i tr Q')/2. The spatial first Ricci matrix is p^2[Q-e(Qe)^t-(Qe)e^t+tr(Q)ee^t]/2. Squaring this complete expression, including its background and volume cross terms, gives the Ricci formula in matching.hessians. Put A=T-2V+UG, B=T-2V+UD, C=V-UD-UG+S. Its coefficients are

    AF p^2 {[(H'+dH^2)E-(A+AD)(H'+(d-2)H^2)]Gamma
      +[(d-4)A/2+dB/2+(d-2)C-AD+AG]H Gamma'
      +[(A+B)/2+C]Gamma''}
    +PF p^4(T-2V+2W+S-UD-UG)Gamma/2.

The background Weyl tensor vanishes. Its second variation is therefore the quadratic linear flat-Weyl form with weight a^(d-3), not the square of only its spatial tracefree part. The raw conformal p^2 density has D''G,D'G',DG'' coefficients

    WA=-2A/(d-1)+4AD/[d(d-1)],
    WB=-4(T-V)+4C/(d-1),
    WC=-2B/(d-1)+4AG/[d(d-1)].

The p^4 coefficient is2(d-2)(T-2V)/(d-1)+2(d-2)W/d-2(d-2)(S-UD-UG)/[d(d-1)]. Integration with the full weight gives the implemented Weyl operator. Finally Riemann^2=Weyl^2+4Ricci^2/(d-1)-2R^2/[d(d-1)].

Setting S=UD=UG=0 recovers every S212 covariant geometric formula; those geometric formulas do not share the amplitude branch error. Independent literal inverse metric, Christoffel, Riemann, contractions and full volume calculations check48 scalar cases at d3,4,5,6 and8 additional noncommuting nonzero-trace cases. The physical Euler Hessian vanishes only after all terms and integrations are retained. Its dimension derivative is not set to zero.
