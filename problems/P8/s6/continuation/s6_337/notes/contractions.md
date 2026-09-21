# Entire local functional before and after subtraction

Define f=Phi, v=dPhi, Y=v.v, H_ab=nabla_a nabla_b Phi,
b=trace H, c=H:H, d=vHv, e=vH^2v and L=b*d-e=L3-L4.
The complete one-pair Wick factors, with I(D) stripped, are

Phi6:                 15 f^4
Phi4Y:                f^4+6 f^2Y
Phi2Y2:               Y^2+(2+4/D) f^2Y
Y3:                   (3+12/D)Y^2
Phi2(L3-L4):          L-2f[(D-2)d+bY]/D+f^2(b^2-c)/D
Y(L3-L4):             (1+4/D)L+Y(b^2-c)/D.

These are OFF-SHELL flat local polynomials before covariantization.
The direct 720-label contraction agrees before using any external
mass shell or momentum-conservation reduction of that polynomial.

One way to derive the fifth row is to note that shifting H by t*g
changes L at first order by (D-2)d+bY. The factor D from the trace
must be retained. Scalar/Hessian contractions give -g I/D,
gradient/gradient contractions give +g I/D, and the Hessian-pair
contractions in L cancel. The last row follows from
Delta_v(Y L)=2(D+4)L+2Y(b^2-c).

Multiply the full factors by I(4-2epsilon) and subtract the pole.
With the overall -1/(16pi^2) retained, each finite factor is
F(4)-2F'(4), giving

Phi6:                 15 f^4
Phi4Y:                f^4+6 f^2Y
Phi2Y2:               Y^2+(7/2) f^2Y
Y3:                   (15/2)Y^2
Phi2(L3-L4):          L-f*d/2-3f*bY/4+3f^2(b^2-c)/8
Y(L3-L4):             (5/2)L+3Y(b^2-c)/8.

Their on-shell four-scalar vertices reproduce ALL SIX frozen S239
finite factors, including its Galileon evanescence. No new flat
subtraction convention is introduced.

The packaged calculation additionally compares the contracted external
Christoffel terms directly to the literal quartic connection of each
F(D). All six identities vanish for generic physical Gram data at
symbolic D, separately from the metric-plus-bubble equality. Thus the
identified quartic representative is not inferred merely from a
four-point flat amplitude. Four nonzero bubble controls and four nonzero
premature-D4 finite errors guard against dropping either contribution.
