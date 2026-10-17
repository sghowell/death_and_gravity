"""Tests of the S4-4 layer driver's bookkeeping (scripts/s44_layer.py):
selection, engine options and the validity check that decides whether a
stored part is skipped.  No region code runs."""
import argparse
import copy
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
import s44_layer as s44  # noqa: E402

CODE_FILE = "problems/P4/FORMULATION.md"


def _args(**kw):
    base = dict(ids=None, regions=None, all_gap=False, beyond_band=False, tau_min=224.0,
                tau_mins=None, norm="layer", combine=False, mv=None,
                piece_ladder_evals=None, piece_final_evals=None, tube_cache=None,
                force=False, max_pieces=0)
    base.update(kw)
    return argparse.Namespace(**base)


@pytest.fixture(scope="module")
def plan():
    return s44.s43.load_plan()


def test_band_membership(plan):
    ps = s44.select(plan, _args())
    ids = [p["id"] for p in ps]
    assert len(ids) == 140
    assert sum(i.startswith("sonic/") for i in ids) == 51
    gap = sorted(p["index"] for p in ps if p["region"] == "gap")
    assert gap == list(range(89))
    # the band is exactly the plan's direct-only pieces
    assert set(ids) == {p["id"] for p in plan["pieces"] if p.get("tq") is False}


def test_beyond_band_selection(plan):
    ps = s44.select(plan, _args(beyond_band=True, regions="tube,deep",
                                ids="tube/full-13[89],tube/full-1[4-9]?,tube/full-2[0-3]?,"
                                    "tube/full-24[0-4],deep/*"))
    assert len(ps) == 236
    tube = sorted(p["index"] for p in ps if p["region"] == "tube")
    assert tube == list(range(138, 245))
    assert not any(p["kind"] == "centre" for p in s44.select(plan, _args(beyond_band=True)))


def test_tau_mins_and_dirs():
    assert s44.tau_mins_of(_args()) == [224.0]
    assert s44.tau_mins_of(_args(tau_mins="224,14")) == [224.0, 14.0]
    with pytest.raises(SystemExit):
        s44.tau_mins_of(_args(tau_mins="14,14"))
    assert s44.layer_dir(224.0).endswith(os.path.join("parts", "layer"))
    assert s44.layer_dir(14.0).endswith(os.path.join("parts", "layer_t14"))


def test_engine_opts_and_validity(plan):
    p = [q for q in plan["pieces"] if q["id"] == "gap/gap-000"][0]
    code = {CODE_FILE: s44.s43.current_hash(CODE_FILE)}
    d = dict(schema=s44.SCHEMA, spec=s44.s43.spec_of(p), code=dict(files=code),
             driver_sha256=s44.DRIVER_SHA256)
    default = s44.engine_opts(_args())
    assert s44.part_valid(d, p, default)                 # no engine record = defaults
    assert not s44.part_valid(dict(d, driver_sha256="8" * 64), p, default)   # other driver
    assert not s44.part_valid(dict(d, code=dict(files=code, write_check=dict(ok=False))),
                              p, default)                # sources changed mid-run
    ext = s44.engine_opts(_args(combine=True, norm="piece", mv="mv3",
                                piece_ladder_evals=200000, piece_final_evals=600000))
    assert not s44.part_valid(d, p, ext)
    d2 = dict(d, engine=ext)
    assert s44.part_valid(d2, p, ext) and not s44.part_valid(d2, p, default)
    d3 = copy.deepcopy(d2)
    d3["spec"]["tq"] = True
    assert not s44.part_valid(d3, p, ext)
    d4 = dict(d2, code=dict(files={CODE_FILE: "0" * 64}))
    assert not s44.part_valid(d4, p, ext)
    kw = s44.engine_kwargs(ext)
    assert kw["norm"] == "piece" and kw["sigma_fn"].__name__ == "sigma_box_mv3"
    assert kw["piece_ladder_evals"] == 200000 and kw["piece_final_evals"] == 600000


def test_work_argv_roundtrip():
    a = _args(regions="gap", tau_mins="224,14", combine=True, norm="piece", mv="mv3",
              piece_ladder_evals=200000, piece_final_evals=600000)
    av = s44.work_argv(a)
    assert av[0] == "work"
    for x in ("--regions=gap", "--tau-mins=224,14", "--norm=piece", "--mv=mv3",
              "--piece-ladder-evals=200000", "--piece-final-evals=600000", "--combine"):
        assert x in av
    assert not any(x.startswith("--tau-min=") for x in av)
    assert "--tau-min=224" in s44.work_argv(_args())


def test_compare_previous_ignores_attempt_times():
    prev = dict(status="pass", certified_bound=-0.5, worst=-0.5, boxes=10, evals=40,
                attempts=[dict(target=-0.75, ok=False, evals=20, time=1.0),
                          dict(target=-0.5, ok=True, evals=20, time=2.0)])
    rec = dict(status="pass", certified_bound=-0.5,
               result=dict(worst=-0.5, boxes=10, evals=40,
                           attempts=[dict(target=-0.75, ok=False, evals=20, time=9.0),
                                     dict(target=-0.5, ok=True, evals=20, time=8.0)]))
    assert s44.compare_previous(prev, rec)["identical"] is True
    rec["result"]["boxes"] = 11
    rc = s44.compare_previous(prev, rec)
    assert rc["identical"] is False and rc["fields"]["boxes"] is False
