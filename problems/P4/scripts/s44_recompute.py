#!/usr/bin/env python3
"""Recompute every layer part on the code with the evaluation-ball fix
(2026-10-17; p4.validated.ballutil, tc_layer_mv's iota-centre): all 546
sonic and tube-kind pieces of the plan, both tau modes, so that every
non-centre piece is carried by a fixed-code layer_t14 part.

Queues (each is one `s44_layer.py work` command; every worker runs them in
the order Q1, Q2, Q4, Q5, Q6, Q3 and the O_EXCL locks of s44_layer.py share
each queue; the locks are per piece, not per tau mode, so Q3, which touches
the pieces of Q2, comes last, long after Q2 has drained):
  Q1 sonic s00-s28           --tau-mins 224,14   (default engine, as before)
  Q2 sonic s29-s50           --tau-min 224       (default engine, as before)
  Q3 sonic s29-s50           --tau-min 14 --norm piece   (as before)
  Q4 gap 000-088             --tau-mins 224,14 --combine --norm piece --mv mv3
                             --piece-ladder-evals 200000 --piece-final-evals 600000
  Q5 tube 138-244, deep      --tau-mins 224,14 --combine --norm piece --mv mv3
                             --piece-ladder-evals 400000 --piece-final-evals 1000000
  Q6 tube 0-137, band, gap 089-106   as Q5
Q1-Q4 use the engine options of the parts they replace; Q5 the options of
s4-4-layer-engine.md section 13; Q6 (pieces without layer parts so far)
the same as Q5.

usage: s44_recompute.py launch --workers 4 [--nice 10] [--tag T]
       s44_recompute.py queues          print the queues
Before launching: s44_layer.py snapshot (old values of every layer part).
Check: s44_layer.py status [--tau-min 14]; s44_layer.py compare.
"""
from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import s43_parts as s43  # noqa: E402

S44 = os.path.join(HERE, "s44_layer.py")
LOGS = os.path.join(s43.P4, "logs", "s44")
MV = ["--combine", "--norm=piece", "--mv=mv3"]
QUEUES = [
    ("Q1", ["--regions=sonic", "--ids=sonic/s0?,sonic/s1?,sonic/s2[0-8]", "--tau-mins=224,14"]),
    ("Q2", ["--regions=sonic", "--ids=sonic/s29,sonic/s3?,sonic/s4?,sonic/s50", "--tau-min=224"]),
    ("Q3", ["--regions=sonic", "--ids=sonic/s29,sonic/s3?,sonic/s4?,sonic/s50", "--tau-min=14",
            "--norm=piece"]),
    ("Q4", ["--regions=gap", "--tau-mins=224,14"] + MV
     + ["--piece-ladder-evals=200000", "--piece-final-evals=600000"]),
    ("Q5", ["--beyond-band", "--regions=tube,deep",
            "--ids=tube/full-13[89],tube/full-1[4-9]?,tube/full-2[0-3]?,tube/full-24[0-4],deep/*",
            "--tau-mins=224,14"] + MV
     + ["--piece-ladder-evals=400000", "--piece-final-evals=1000000"]),
    ("Q6", ["--beyond-band", "--regions=tube,band,gap",
            "--ids=tube/full-0??,tube/full-1[0-2]?,tube/full-13[0-7],band/*,gap/gap-089,"
            "gap/gap-09?,gap/gap-10?",
            "--tau-mins=224,14"] + MV
     + ["--piece-ladder-evals=400000", "--piece-final-evals=1000000"]),
]


ORDER = ("Q1", "Q2", "Q4", "Q5", "Q6", "Q3")


def chain_cmd(log):
    parts = []
    qd = dict(QUEUES)
    for name in ORDER:
        av = qd[name]
        cmd = " ".join([sys.executable, S44, "work"] + ["'%s'" % a for a in av])
        parts.append("echo '== %s' >> %s; %s >> %s 2>&1" % (name, log, cmd, log))
    return " ; ".join(parts)


def cmd_launch(args):
    tag = args.tag or time.strftime("s44-recompute-%Y%m%d-%H%M%S")
    os.makedirs(LOGS, exist_ok=True)
    workers = []
    for k in range(args.workers):
        log = os.path.join(LOGS, "%s-w%d.log" % (tag, k))

        def pre():
            os.nice(args.nice)
            signal.signal(signal.SIGHUP, signal.SIG_IGN)
        p = subprocess.Popen(["nohup", "sh", "-c", chain_cmd(log)], cwd=s43.ROOT,
                             stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                             stderr=subprocess.DEVNULL, start_new_session=True,
                             preexec_fn=pre, env=dict(os.environ, PYTHONUNBUFFERED="1"))
        workers.append(dict(pid=p.pid, log=s43.rel(log)))
    rec = dict(tag=tag, host=os.uname().nodename, started=s43.now(), nice=args.nice,
               workers=workers, order=list(ORDER),
               queues=[dict(name=n, argv=["work"] + av) for n, av in QUEUES],
               git_head=s43.git_head()[0],
               s44_layer_sha256=s43.sha256_file(S44),
               note="each worker: sh -c running the queues in order (pid = the sh)")
    s43.write_json_atomic(os.path.join(LOGS, "launch-%s.json" % tag), rec)
    for w in workers:
        print("pid %d  log %s" % (w["pid"], w["log"]))


def cmd_queues(args):
    qd = dict(QUEUES)
    for name in ORDER:
        print(name, " ".join(qd[name]))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    q = sub.add_parser("launch")
    q.add_argument("--workers", type=int, required=True)
    q.add_argument("--nice", type=int, default=10)
    q.add_argument("--tag")
    q.set_defaults(fn=cmd_launch)
    q = sub.add_parser("queues")
    q.set_defaults(fn=cmd_queues)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
