# Counterterm structure at arbitrary graph loop order

Consider connected graphs with V_G cubic H Phi^2 vertices,
V_4 quartic Phi^4 vertices, internal line counts I_phi,I_H
and external counts E_phi,E_H. Exact incidences and topology are

    2 I_phi+E_phi=2 V_G+4 V_4,
    2 I_H+E_H=V_G,
    L=I_phi+I_H-V_G-V_4+1.

The superficial momentum degree in four dimensions is

    d=4L-2I_phi-2I_H=4-E_phi-E_H-V_G.

Every graph and subgraph obeys E_phi even and V_G>=E_H,
with V_G-E_H even. If d>=0 then E_phi+E_H+V_G<=4.
This proves that the finite exact enumeration in powercount.py
exhausts the potentially divergent structures for UNBOUNDED
loop order and arbitrary numbers of quartic vertices.

For H-dependent terms only (E_phi,E_H)=(0,1),(0,2),(2,1)
remain, with maximum momentum degrees 2,0,0 respectively.
The one-point external momentum is zero by translation
invariance; a derivative on a single H is a boundary term.
Thus the needed structures are the H tadpole, H^2 mass and
H Phi^2 vertex. There is no divergent H kinetic term,
H^3, H^4 or H^2 Phi^2. Constant finite H kinetic
renormalization is still allowed and remains Gaussian.

Pure-light terms have only the usual mass, kinetic and
quartic structures, in addition to a vacuum constant.
Using standard local perturbative subtraction of subdivergences,
the same counting applies inductively to subgraphs. Mass or
tadpole insertions can only lower the degree; renormalizable
kinetic and quartic counterterm insertions do not enlarge it.
Some listed diagrams may vanish; that cannot add a new structure.

This proves a power-counting counterterm-basis closure. It
does not compute finite coefficients, prove convergence of
the perturbation series, remove a regulator nonperturbatively
or establish the complex-energy hypotheses of V. In particular,
a Gaussian H path integral leaves interacting phi loops and
does not amount to a free quantum theory.
