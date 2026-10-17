#!/usr/bin/env python3
"""S4-6 review: the default merge and `merge --substitute-layer-t14` are
unchanged by the S4-6 edit of s43_parts.py.

Runs cmd_merge of main's s43_parts.py (`git show main:...`, executed
in-process with __file__ set to the real script path, so its paths resolve
to this checkout) and of this branch's, on the same committed parts, in both
modes, with the outputs redirected to a temporary directory (no tracked file
written).  Reports every JSON path that differs between main and the branch.
The S4-6 edit is expected to change only the merge-output wording
(scope.full_strip_tau and scope.xi_open_max replaced; scope.tau_target,
tau_by_carrier, open_at_centre_end added), add centre_substitution, the
phase_settings note and provenance fields; never status, the cover, detail,
missing, failed, substitution or strip_exclusion_certified.

usage (repository root): PYTHONPATH=problems/P4/src python .../merge_unchanged.py [out.json]
"""
import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
SCRIPT = os.path.join(ROOT, "problems", "P4", "scripts", "s43_parts.py")


def load(src, name):
    mod = types.ModuleType(name)
    mod.__file__ = SCRIPT
    sys.modules[name] = mod
    exec(compile(src, SCRIPT, "exec"), mod.__dict__)
    return mod


def run(mod, tmp, layer, tag):
    mod.CERT = os.path.join(tmp, tag + "_complete.json")
    mod.CERT_PARTIAL = os.path.join(tmp, tag + "_partial.json")
    args = types.SimpleNamespace(accept_code=None, accept_reason=None,
                                 substitute_layer_t14=layer, substitute_centre=False)
    with contextlib.redirect_stdout(io.StringIO()):
        mod.cmd_merge(args)
    p = mod.CERT if os.path.exists(mod.CERT) else mod.CERT_PARTIAL
    return json.load(open(p))


def walk(a, b, path=""):
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b), key=str):
            out += walk(a.get(k, "<absent>"), b.get(k, "<absent>"), path + "." + str(k))
        return out
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        out = []
        for i, (x, y) in enumerate(zip(a, b)):
            out += walk(x, y, path + "[%d]" % i)
        return out
    return [] if a == b else [path]


def main(out):
    src_main = subprocess.run(["git", "-C", ROOT, "show", "main:problems/P4/scripts/s43_parts.py"],
                              capture_output=True, text=True, check=True).stdout
    src_head = open(SCRIPT).read()
    tmp = tempfile.mkdtemp(prefix="s46_unchanged_")
    m_main = load(src_main, "s43_main")
    res = {}
    for layer in (False, True):
        res[("main", layer)] = run(m_main, tmp, layer, "main_%d" % layer)
    m_head = load(src_head, "s43_head")
    for layer in (False, True):
        res[("head", layer)] = run(m_head, tmp, layer, "head_%d" % layer)
    rec = dict(schema="review-s46-merge-unchanged/1", modes={})
    for layer in (False, True):
        a, b = res[("main", layer)], res[("head", layer)]
        d = walk(a, b)
        tops = sorted({p.split("[")[0].split(".")[1] for p in d})
        scope = sorted({".".join(p.split(".")[:3]) for p in d if p.startswith(".scope")})
        key = dict(status=(a["status"], b["status"]),
                   certified=(a["strip_exclusion_certified"], b["strip_exclusion_certified"]),
                   missing=(len(a["missing"]), len(b["missing"])),
                   substituted=(a["substitution"]["count"], b["substitution"]["count"]),
                   gaps=(a["scope"]["tau_gt_tstar_uncovered_x"],
                         b["scope"]["tau_gt_tstar_uncovered_x"]))
        same_core = all(walk(a.get(k), b.get(k)) == [] for k in
                        ("ok", "status", "cover", "detail", "missing", "failed", "substitution",
                         "strip_exclusion_certified", "checks", "checks_valid_for_code",
                         "code_conflicts", "pieces", "worst_margin", "layer",
                         "layer_join_rejections", "open_pieces"))
        rec["modes"]["substitute_layer_t14" if layer else "default"] = dict(
            differing_top_level_keys=tops, differing_scope_keys=scope,
            core_identical=same_core, key=key, n_paths=len(d))
    rec["ok"] = all(v["core_identical"] for v in rec["modes"].values())
    json.dump(rec, open(out, "w"), indent=1)
    print(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "merge_unchanged.json"))
