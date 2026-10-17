#!/usr/bin/env python3
"""S4-6 review: provenance of the centre parts, the identity artifact and the
pinned code, against git (no p4 code imported).

  * every recorded source hash of the 60 centre parts, the identity artifact
    and the 1092 layer_t14 / layer parts equals the file at HEAD;
  * the 26 sources the layer parts pin and scripts/s44_layer.py are
    byte-identical on main and at HEAD;
  * problems/P4/src is unchanged between the commits the parts and the
    identity record (e9cc5086, da7a7f28) and HEAD; s46_centre.py too;
  * the centre parts' provenance fields (git head, clean tree, no late
    module, passing write check, driver) and the sha256 / aggregates the
    notes quote.

usage (repository root): python problems/P4/results/theorem_c/review_s46/provenance_check.py [out.json]
"""
import glob
import hashlib
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
P4 = os.path.join(ROOT, "problems", "P4")
TC = os.path.join(P4, "results", "theorem_c")


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def git(*a):
    return subprocess.run(["git", "-C", ROOT] + list(a), capture_output=True, text=True,
                          check=True).stdout


def git_sha(rev, rel):
    out = subprocess.run(["git", "-C", ROOT, "show", "%s:%s" % (rev, rel)],
                         capture_output=True, check=True).stdout
    return hashlib.sha256(out).hexdigest()


def aggregate(d):
    files = sorted(glob.glob(os.path.join(d, "**", "*.json"), recursive=True),
                   key=lambda p: os.path.relpath(p, ROOT).encode())
    lines = "".join("%s  %s\n" % (sha(p), os.path.relpath(p, ROOT)) for p in files)
    return hashlib.sha256(lines.encode()).hexdigest(), len(files)


def main(out):
    rec = dict(schema="review-s46-provenance/1", head=git("rev-parse", "HEAD").strip())
    cur = {}

    def current(rel):
        if rel not in cur:
            p = os.path.join(ROOT, rel)
            cur[rel] = sha(p) if os.path.exists(p) else None
        return cur[rel]

    s46 = os.path.join(P4, "scripts", "s46_centre.py")
    s44 = os.path.join(P4, "scripts", "s44_layer.py")
    # centre parts
    bad, heads, drivers = [], set(), set()
    for t in (14, 224):
        for p in sorted(glob.glob(os.path.join(TC, "parts", "centre_t%d" % t, "centre", "*.json"))):
            d = json.load(open(p))
            c = d["code"]
            stale = [r for r, h in c["files"].items() if current(r) != h]
            heads.add(c["git_head"])
            drivers.add(d["driver_sha256"])
            if stale or c["git_dirty_src"] is not False or c["late_modules"] \
                    or not c["write_check"]["ok"] or d["driver_sha256"] != sha(s46) \
                    or d["status"] != "pass":
                bad.append(dict(part=os.path.relpath(p, ROOT), stale=stale))
    rec["centre_parts"] = dict(bad=bad, git_heads=sorted(heads), drivers=sorted(drivers),
                               s46_current=sha(s46))
    # identity artifact
    idp = os.path.join(TC, "centre_nb_identity.json")
    idr = json.load(open(idp))
    rec["identity"] = dict(sha256=sha(idp), git_head=idr["code"]["git_head"],
                           dirty=idr["code"]["git_dirty_src"],
                           stale=[r for r, h in idr["code"]["files"].items() if current(r) != h],
                           driver_current=idr["driver_sha256"] == sha(s46), ok=idr["ok"])
    # layer parts
    lbad, pinned, ldrv = [], set(), set()
    nlayer = 0
    for dname in ("layer_t14", "layer"):
        for p in glob.glob(os.path.join(TC, "parts", dname, "**", "*.json"), recursive=True):
            d = json.load(open(p))
            if "code" not in d:
                continue
            nlayer += 1
            files = d["code"]["files"]
            pinned |= set(files)
            ldrv.add(d.get("driver_sha256"))
            stale = [r for r, h in files.items() if current(r) != h]
            if stale:
                lbad.append(dict(part=os.path.relpath(p, ROOT), stale=stale))
    pin = sorted(pinned)
    same_main = {r: git_sha("main", r) == current(r) for r in pin}
    rec["layer_parts"] = dict(count=nlayer, stale=lbad[:10], n_stale=len(lbad),
                              pinned_sources=len(pin),
                              pinned_identical_on_main=all(same_main.values()),
                              pinned_differing=[r for r, v in same_main.items() if not v],
                              drivers=sorted(x for x in ldrv if x),
                              s44_current=sha(s44),
                              s44_identical_on_main=git_sha(
                                  "main", "problems/P4/scripts/s44_layer.py") == sha(s44))
    # code history
    rec["src_diff_main_HEAD"] = git("diff", "--stat", "main", "HEAD", "--", "problems/P4/src")
    rec["src_diff_e9cc5086_HEAD"] = git("diff", "--stat", "e9cc5086", "HEAD", "--",
                                        "problems/P4/src", "problems/P4/scripts/s46_centre.py")
    rec["src_diff_da7a7f28_HEAD"] = git("diff", "--stat", "da7a7f28", "HEAD", "--",
                                        "problems/P4/src", "problems/P4/scripts/s46_centre.py")
    rec["src_diff_fda58897_HEAD"] = git("diff", "--stat", "fda58897", "HEAD", "--",
                                        "problems/P4/src")
    # quoted hashes
    q = dict(
        centre_t14=aggregate(os.path.join(TC, "parts", "centre_t14")),
        centre_t224=aggregate(os.path.join(TC, "parts", "centre_t224")),
        layer_t14=aggregate(os.path.join(TC, "parts", "layer_t14")),
        identity=sha(idp), tc_commutator=sha(os.path.join(TC, "tc_commutator.json")),
        sampling=sha(os.path.join(TC, "centre_sampling.json")),
        tc_layer_centre=sha(os.path.join(P4, "src", "p4", "validated", "tc_layer_centre.py")),
        s46_centre=sha(s46))
    want = dict(centre_t14="4ac8a7b00f9904d56dd5a26e7b656fd6387f1457614be7467dc6cad61e1f17d7",
                centre_t224="ba424c3f78ea8d4a001072dd07025b7c6cd56acbfce12c5304ad14cf7f9fc836",
                layer_t14="e9835a41fb2dbcc718e6beff9d6aa89b3b50b7874c996baee10fa173d0f5ff6c",
                identity="0fad868ba62fdecb836e6f8a78de23bd3a87f096e50bd117b3cb0475c29bd505",
                tc_commutator="6c228a73abb42c2a2a17ef0c59f4e6565e5c4381b19b7d1dafc15556150fdbab",
                sampling="a740dd19abea81c74193caf7f3f788cad11f57f1597c4ddbe476c0382f789c06",
                tc_layer_centre="ec32ccd6", s46_centre="735336cc")
    rec["quoted"] = {k: dict(computed=(v if isinstance(v, str) else v[0]),
                             files=(None if isinstance(v, str) else v[1]), quoted=want[k],
                             match=(v if isinstance(v, str) else v[0]).startswith(want[k]))
                     for k, v in q.items()}
    rec["ok"] = (not bad and not rec["identity"]["stale"] and rec["identity"]["driver_current"]
                 and not lbad and rec["layer_parts"]["pinned_identical_on_main"]
                 and rec["layer_parts"]["s44_identical_on_main"]
                 and not rec["src_diff_e9cc5086_HEAD"].strip()
                 and not rec["src_diff_da7a7f28_HEAD"].strip()
                 and all(v["match"] for v in rec["quoted"].values()))
    json.dump(rec, open(out, "w"), indent=1)
    print(json.dumps({k: rec[k] for k in ("head", "ok")}, indent=1))
    print("centre parts bad:", len(bad), "heads", rec["centre_parts"]["git_heads"])
    print("identity:", rec["identity"])
    print("layer parts:", nlayer, "stale", len(lbad), "pinned", len(pin),
          "identical on main", rec["layer_parts"]["pinned_identical_on_main"],
          "s44 identical", rec["layer_parts"]["s44_identical_on_main"])
    print("src diff main..HEAD:\n" + rec["src_diff_main_HEAD"])
    print("src diff e9cc5086..HEAD: %r" % rec["src_diff_e9cc5086_HEAD"])
    print("src diff da7a7f28..HEAD: %r" % rec["src_diff_da7a7f28_HEAD"])
    print("src diff fda58897..HEAD:\n" + rec["src_diff_fda58897_HEAD"])
    for k, v in rec["quoted"].items():
        print("  %-16s match=%s %s (%s files)" % (k, v["match"], v["computed"][:16], v["files"]))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "provenance_check.json"))
