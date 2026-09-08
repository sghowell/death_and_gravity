"""Independent temporal solve, canonical force and normalization controls."""
import sympy as sp
from p8_aligned_response import forcing as f


def test_actual_forced_equation_and_readout_identities():
    assert all(value == 0 for value in f.checks().values())


def test_full_temporal_solve_independent_of_stated_solution():
    data = f.temporal()
    solved = sp.solve(sp.diff(data["L"], data["t"]), data["t"])
    assert len(solved) == 1
    assert sp.factor(solved[0]-data["solution"]) == 0
    assert sp.factor(data["L"].subs(data["t"], solved[0])-data["reduced"]) == 0
    assert solved[0].subs(data["q"], 0) == data["S"]


def test_unweighted_source_and_omitted_normalization_are_wrong():
    data = f.canonical()
    g, source = data["g"], data["source"]
    assert sp.expand(data["force"]-source) != 0
    assert sp.expand(data["force"]-g*sp.diff(source, f.u)) == 2*sp.diff(g, f.u)*source


def test_preparation_requires_source_jets_not_merely_source_value():
    data = f.canonical()
    source = data["source"]
    S0, S1, S2 = sp.symbols("S0 S1 S2", real=True)
    jets = {source: S0, sp.diff(source, f.u): S1, sp.diff(source, f.u, 2): S2}
    assert data["force"].xreplace(jets).subs(S0, 0) == data["g"]*S1
    assert data["force_prime"].xreplace(jets).subs({S0: 0, S1: 0}) == data["g"]*S2


def test_nonzero_second_order_source_is_imported_not_a_spectator():
    from p8_affine_aligned import alignment
    assert alignment.rolling()["second_order"] != 0
    assert alignment.quadratic_mass()["third_order_source_coupling"] != 0
