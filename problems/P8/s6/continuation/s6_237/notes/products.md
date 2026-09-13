# Complete product remainder after the fourth amplitude order

Let w=1/n<=1/64. The exact vertex remainder is

r_a=w5[s4/(1-sw)+16/(1-2w)²].

Because |s|<3,

|r_a|/w5<=81/(1-3/64)+16/(1-2/64)²
             =5981248/58621<128.

For the four retained vertex coefficients use the absolute bounds[2,7,13,27]. The polynomial norm divided by w is at most553819/262144. Adding128/64^4 still leaves the full vertex norm below3w.

The four retained triangle coefficients obey |c_j|<=2M*32^(j-1), j=1,...,4, by the complete P,Q Cauchy coefficient bounds. Therefore the full retained triangle polynomial has modulus below4Mw.

## Bubble

In the square of the four-term vertex polynomial, sum EVERY product with total order at least6. Its full positive coefficient majorant, including its later powers of w, is

(1/2) sum_{i,j=1..4,i+j>=6} a_i a_j (1/64)^(i+j-6)
 =2286169/8192<512.

The cross term with the exact r_a is below384w6, and r_a²/2 is below w6. With |B|<2, the entire bubble remainder is below2048w6, hence also below2048Mw6. Terms of degree above6 have not been omitted.

## Triangle

Write the full product2 a Cbar as its four-order polynomial part plus three complete remainder groups.

The exact full vertex times R_C is bounded by6*2^22 M w6. The exact r_a times the retained triangle polynomial is bounded by1024 M w6.

For the product of the two retained polynomials, enumerate every term of total degree at least6. Its full norm coefficient is

2 sum_{i,j=1..4,i+j>=6} a_i[2*32^(j-1)](1/64)^(i+j-6)
 =1003424<2^20,

multiplying M w6. This sum includes the w7 and w8 cross products. These groups partition the full error without dropping or double-counting r_a R_C: that product is already in the exact full vertex times R_C.

## Two ordered boxes and all three channels

Each box remainder is below3*2^23 M w6, and there are two per channel. Combining ALL terms gives a per-channel majorant

[2048+6*2^22+1024+2^20+2*3*2^23]M w6<2^27 M w6.

All three channels have constant3*2^27*500<10^12. After restoring the common normalization,

abs(R6_amplitude)<10^12 g4/(16pi² n6)

uniformly on the complete closed complex bidisk. This estimate concerns the full first-loop amplitude minus its explicit n^-2 through n^-5 terms. It is not an omitted-loop estimate.
