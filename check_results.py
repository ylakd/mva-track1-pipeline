#!/usr/bin/env python3
"""Compare a fresh pipeline run with the numbers reported in the methods report.

    python check_results.py [--workdir DIR]

A difference is not necessarily an error (live services and the 'latest' HPO release can change),
but each one must be explained before the report is submitted.
"""
import argparse, json, os, sys
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--workdir", default=None)
    args = ap.parse_args()
    if args.workdir:
        os.environ["MVA_WORKDIR"] = args.workdir
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from mvapipe import config as cfg
    with open(Path(__file__).with_name("expected_results.json")) as fh:
        exp = json.load(fh)
    wd, bad = cfg.workdir(), 0

    def load(name):
        p = wd / name
        if not p.exists():
            return None
        with open(p) as fh:
            return json.load(fh)

    print(f"{'item':58} {'expected':>12} {'observed':>12}  status")
    for fname, items in exp.items():
        if fname.startswith("_") or fname in ("text", "sets"):
            continue
        data = load(fname)
        for key, (val, tol) in items.items():
            obs = None if data is None else data.get(key)
            ok = obs is not None and abs(obs - val) <= tol
            bad += not ok
            print(f"{fname + ':' + key:58} {val:>12} {str(obs)[:12]:>12}  {'PASS' if ok else 'DIFF'}")
    for fname, items in exp.get("text", {}).items():
        data = load(fname) or {}
        for key, val in items.items():
            ok = data.get(key) == val
            bad += not ok
            print(f"{fname + ':' + key:58} {val:>12} {str(data.get(key)):>12}  {'PASS' if ok else 'DIFF'}")
    for fname, items in exp.get("sets", {}).items():
        data = load(fname) or {}
        for key, val in items.items():
            ok = sorted(data.get(key, [])) == sorted(val)
            bad += not ok
            print(f"{fname + ':' + key:58} {'set':>12} {'match' if ok else 'differs':>12}  {'PASS' if ok else 'DIFF'}")
    print("\nAll items match." if not bad else f"\n{bad} item(s) differ: explain each before submitting.")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
