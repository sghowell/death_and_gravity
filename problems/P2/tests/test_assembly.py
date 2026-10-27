"""S5 assembly (P2-4 at P_B; notes/s5-plan.md section 3, notes/s5-assembly.md).

Fast (~5 s): the bookkeeping check on the committed certificates passes; every negative control fails on
its named item; the stored certificates/assembly/assembly_PB.json matches a recomputation and pins the
current inputs; the ingredient levels used for P2-4's level match CLAIMS.md.
"""
import json

import pytest

from p2 import assembly


@pytest.fixture(scope="module")
def run():
    # round trip through JSON so that the comparison with the stored file sees the same key types
    return json.loads(json.dumps(assembly.assemble(write=False), ensure_ascii=False))


def test_check_passes(run):
    r = run["result"]
    assert r["pass"] and r["failed"] == []
    s = r["summary"]
    assert s["d_delta_over_kappa_minus"] == "11/20"
    assert s["l0"] == 20
    per_l = s["l_le_19"]["per_l"]
    assert sorted(int(k) for k in per_l) == list(range(20))
    assert all(per_l[str(l)]["winding"] == (1 if l == 0 else 0) for l in range(20))
    assert s["l_le_19"]["b"] == ["6280879/268435456"] and s["l_le_19"]["e"] == ["4567911/1073741824"]
    assert s["levels"]["minimum"] == "WRITTEN_PROOF"
    assert s["hv17"]["exponent"] == "21/20"
    assert s["cf17"]["two_min_over_kappa_minus"] == "11/10"
    for name in ("C3.cover_replay_identical", "C3.no_l_gap", "C4.E_plus_below_strip", "C4.E_c_below_strip",
                 "C5.D4_derivative_nonzero", "C6.beta_upper", "C7.alpha_lt_gamma0", "C8.km_lt_2kp"):
        assert r["checks"][name]["ok"], name


def test_negative_controls_fail(run):
    ctr = run["negative_controls"]
    assert [c["id"] for c in ctr] == ["N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8a", "N8b"]
    for c in ctr:
        assert c["fails_as_expected"], c["id"]
        assert set(c["expected_failures"]) <= set(c["failed"]), c["id"]
    # N1 is the tight one: b just below d_delta fails on the strip-depth check alone
    n1 = ctr[0]
    assert "C2.l0.b_ge_d" in n1["failed"] and len(n1["failed"]) <= 2


def test_stored_certificate(run):
    assert assembly.OUT.exists(), "run python -m p2.assembly"
    stored = json.loads(assembly.OUT.read_text())
    assert stored["pass"] is True
    assert stored["inputs_sha256"] == run["inputs_sha256"], "an input certificate changed since the assembly ran"
    assert stored["result"] == run["result"]
    assert stored["negative_controls"] == run["negative_controls"]


def test_levels_match_claims():
    claims = assembly.ROOT.parents[1] / "CLAIMS.md"
    if not claims.exists():
        pytest.skip("CLAIMS.md not found")
    rows = {}
    for line in claims.read_text().splitlines():
        if line.startswith("| P2-"):
            cols = [c.strip() for c in line.split("|")]
            rows[cols[1]] = cols[4]
    for key, (level, _) in assembly.INGREDIENTS.items():
        row = rows[key.replace("'", "′")]
        assert row.split()[0].rstrip(".") == level, key
