# S4-6 review checks (2026-10-17)

The scripts and outputs of the independent review
`docs/review-2026-10-17-p4-s4-6-centre-adversarial.md` (its section 9 lists what
each one checks and what it found). They ran in the review worktree, a git
checkout of branch `p4/theorem-c-centre-review` made from `p4/theorem-c-centre`
at `36c0c802`, with the gitignored tube caches copied into
`problems/P4/data/cache/` (content hashes equal to the cache records').
Python 3.11.14, python-flint 0.9.0, sympy 1.14.0, at nice 10.

From the repository root, with `PY=.venv/bin/python`, `PP=PYTHONPATH=problems/P4/src`
and `R=problems/P4/results/theorem_c/review_s46`:

| script | command | output | what |
|---|---|---|---|
| `identity_independent.py` | `$PP $PY $R/identity_independent.py` | `identity_independent.json` | U6 built as U1/Dp^6 from the unreduced pair and equal to build_pairs' formula; symbolic series substitution, odd coefficients included: the y^8 and y^9 coefficients of nB's kappa^1..7 coefficients vanish for every background with A = 1 + O(y^2); three controls |
| `raise_inputs_check.py` | `$PP $PY $R/raise_inputs_check.py c00,c01,c15,c28,c29` | `raise_inputs_check.json` | the certified (a*, mu*) balls; odd series coefficients exactly 0; A = 1 + y^2 Ahat; per piece, the dropped TM coefficients lie in the exact range (ve >= 25) and contain 0 in all three jet slots; raise_nb is a drop of exactly two; the (tt dT)' raise never fires |
| `sympy_forms.py` | `$PY $R/sympy_forms.py` | `sympy_forms.json` | T(lambda) = lambda^2 S/tau^2 for a free lambda (control nonzero); the closed forms of both regimes; the Neumann remainder and the fold |
| `gv_stress.py` | `$PP GV_HQMAX=-0.3 $PY $R/gv_stress.py 20000` | `gv_stress.json` | randomized containment of GV value, product, power and inverse; negative exponents |
| `cover_replay.py` | `$PP $PY $R/cover_replay.py c00:14 c00:224 c29:224 c29:14 c27:14 c28:14` | `trees/*.json.gz`, `logs/cover_replay_*.txt` | the parts' covers replayed (tc_layer.certify_piece, copied with logging) with the box tree logged; counts bit-identical to the committed parts |
| `cover_verify.py` | `$PY $R/cover_verify.py --controls $R/trees/c29_t14.json.gz $R/trees/*.json.gz --out $R/cover_verify.json` | `cover_verify.json` | exact (Fraction) proof, from each logged tree, that the accepted boxes cover piece x [0, 6] x [tau_min, inf); five corruption controls |
| `direct_chain.py` | (module) | — | an independent pointwise S from the assembly's definitions, at the midpoint background of the certified balls |
| `point_check.py` | `$PP $PY $R/point_check.py $R/trees/c00_t14.json.gz $R/trees/c00_t224.json.gz $R/trees/c29_t224.json.gz $R/trees/c29_t14.json.gz $R/trees/c27_t14.json.gz $R/trees/c28_t14.json.gz --grid 300 --out $R/point_check.json` | `point_check.json`, `logs/point_check.txt` | direct values at points inside every accepted box lie in the box's enclosure |
| `adversarial_points.py` | `$PP $PY $R/adversarial_points.py` | `adversarial_points.json` | extreme and boundary points (y down to 1e-40, tau up to 1e40, grid edges), checked in every accepted box that contains them |
| `cover_x_full.py` | `$PY $R/cover_x_full.py` | `cover_x_full.json` | independent exact x-cover of all 576 pieces; the junction at the deep end |
| `merge_tamper.py` | `$PP $PY $R/merge_tamper.py` | `merge_tamper.json` | 50 end-to-end merge scenarios (tampered parts, identity artifact, checks.json), in-process; no tracked file is written |
| `merge_unchanged.py` | `$PP $PY $R/merge_unchanged.py` | `merge_unchanged.json` | main's s43_parts.py against the branch's: the default merge and --substitute-layer-t14 have identical core fields |
| `provenance_check.py` | `$PY $R/provenance_check.py` | `provenance_check.json` | recorded hashes current; pinned sources identical to main; the quoted hashes and aggregates |
| `notes_check.py` | `$PY $R/notes_check.py` | `notes_check.json` | the numbers of notes/s4-6-centre.md section 5 against the data; the assembly edits made in place |
| `recompute_compare.py` | after `s46_centre.py work --force` (below): `$PY $R/recompute_compare.py` | `recompute_compare.json`, `merge_recomputed.json` | the 60 recomputed parts against the committed ones; the merge on them |
| `s46_centre.py identity` | `$PP $PY problems/P4/scripts/s46_centre.py identity --out $R/identity_replay.json` | `identity_replay.json` | the identity artifact regenerated: identical in content to the committed one |

**The recompute.** It ran as three disjoint detached groups:

    $PY problems/P4/scripts/s46_centre.py work --pieces c29,c00,c01,c02 --tau-mins 224,14 --max-evals 150000 --force
    $PY problems/P4/scripts/s46_centre.py work --pieces c28,c27,c26,c25,c24,c03,c04,c05,c06,c07,c08,c09 --tau-mins 224,14 --max-evals 150000 --force
    $PY problems/P4/scripts/s46_centre.py work --pieces c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23 --tau-mins 224,14 --max-evals 150000 --force

Their stdout is in `logs/recompute_g1.txt`, `g2.txt` and `g3.txt`. They
overwrote the tracked parts, and the merge in `recompute_compare.py`
overwrote `tc_commutator.json`. After the comparison all of these were
restored with

    git checkout -- problems/P4/results/theorem_c/parts/centre_t14 problems/P4/results/theorem_c/parts/centre_t224 problems/P4/results/theorem_c/tc_commutator.json

**`logs/`** holds captured stdout as `.txt`, because the repository ignores
`*.log`:
- `pytest_p4_full.txt`: the P4 test suite, with the known-slow
  `test_tc_commutator.py::test_centre_first_piece` deselected;
- `pytest_p4_full_attempt1_stopped.txt`: the first attempt, stopped in that
  test;
- the recompute workers, the replays, and the point checks.
