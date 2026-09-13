# Full diagram budget and the actual uniform matching window

Write w=1/n and k=s-4. For S>=4 and S/n<=1/8, n>=32. The exact fixed-vertex remainder from the previous note obeys

|r_a|<=w³[S²/(1-Sw)+4(1+4w)/(1-2w)²]
     <=4S²w³.

The constant estimate follows from1/(1-Sw)<=8/7, w<=1/32, S²>=16 and8/7+8/25<4. Consequently

|A_s/g²|<=2w+2Sw²+4S²w³<3w,
|k|<=2S.

Let M=1000 be the simultaneous bound on |B| and integral |H0|+3 from the preceding note. This is a maximum of two bounds; it does not assume |B| is bounded by integral |H0| itself.

The COMPLETE bubble error after its n^-2 and n^-3 terms is at most16 M S²w4. Indeed its three remainder pieces have constants2,12 and8(Sw)²<=1/8, whose sum113/8 is below16.

For the full triangle2(A_s/g²)Cbar, group the remaining terms as the exact vertex times R_C, r_a times the two displayed Cbar terms, and the cross product of their linear jets. Their respective conservative constants are384,12 and14, all multiplying M S²w4. Thus the entire triangle error is below512 M S²w4. This retains all terms containing r_a; it does not apply the vertex truncation to the full loop.

Each complete ordered box remainder is at most256 M S²w4. There are TWO per channel. Therefore the channel sum is bounded by

(16+512+2*256)M S²w4<2048 M S²w4.

All three channels together have constant3*2048*1000<10^7. Reinstating g4/(16pi²) proves

|R_base|<10^7 g4 S²/(16pi²n4).

This is a bound on the entire complex coefficient, including all real and imaginary contributions, not solely on the first elastic cut.

## Same finite value subtraction, not another parameter adjustment

The symmetric point has the same universal part and obeys the preceding bound with S=4. The existing OS4 contact therefore gives

|A1_OS4|<10^7 g4(S²+16)/(16pi²n4)
       <=2*10^7 g4 S²/(16pi²n4).

For the original full massive tree, S233 proves

A_original>=lambda(S-2)²>=lambda S²/4>0.

The last inequality follows from4(S-2)²-S²=(S-4)(3S-4)>=0. Using lambda=g²/(2D³), n=D+2 and pi²>9,

|A1_OS4|/A_original
 <10^7 g²D³/(pi²n4)
 <10^7 g²/(9n)
 <10^-199.

Every last parameter inequality is checked as an exact rational comparison. In particular S_max=10^196<n/8 for the unchanged g²=1/2^26 and D=10^200/512.

For the full energy window2<=E<=10^98 and every physical angle, S233's tree error is below1/60. The triangle inequality gives

|A_tree,V2S+A1_OS4-A_original|/A_original
 <1/60+10^-199<1/59.

This last expression is explicitly the tree-plus-first-loop truncation with hbar=1 as a bookkeeping convention. The bound on its first-loop term does not estimate the omitted O(hbar²) part. The positive classical contact margin and first-order external normalization remain the ones already proved in S235/S234, not new all-order vacuum statements.
