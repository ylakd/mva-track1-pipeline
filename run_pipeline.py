#!/usr/bin/env python3
"""Run the Track 1 pipeline end to end, or selected stages.

    python run_pipeline.py                      # all stages 1-12
    python run_pipeline.py --stages 2-6         # a range
    python run_pipeline.py --stages 8,9,10,11   # supplementary analyses only
    python run_pipeline.py --workdir /tmp/mva

Stages: 1 data+provenance | 2 QC | 3 coding | 4 annotation | 5 recessive scan | 6 phenotype ranking +
permutations | 7 verification | 8 [supp] locus review | 9 [supp] NMD | 10 [supp] gnomAD |
11 [supp] allele-fraction screen | 12 submission CSV | 13 [test] sensitivity | 14 [test] benchmark
"""
from __future__ import annotations

import argparse, json, os, platform, sys, time


def parse_stages(spec):
    if spec == "all":
        return list(range(1, 15))
    out = []
    for part in spec.split(","):
        if "-" in part:
            a, b = part.split("-"); out += list(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return sorted(set(out))


def versions():
    v = {"python": platform.python_version()}
    for m in ("pandas", "numpy", "cyvcf2", "pyarrow", "requests", "huggingface_hub", "docx"):
        try:
            v[m] = getattr(__import__(m), "__version__", "n/a")
        except Exception:
            v[m] = "not installed"
    return v


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stages", default="all")
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)
    if args.workdir:
        os.environ["MVA_WORKDIR"] = args.workdir
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from mvapipe import config as cfg
    from mvapipe import stage01_03_data_qc_coding as s1
    from mvapipe import stage04_05_annotate_scan as s4
    from mvapipe import stage06_phenotype as s6
    from mvapipe import stage07_verify as s7
    from mvapipe import stage08_locus as s8
    from mvapipe import stage09_10_nmd_gnomad as s9
    from mvapipe import stage11_baf as s11
    from mvapipe import stage12_submission as s12
    from mvapipe import stage13_14_robustness as s13

    registry = {
        1: ("data + provenance", lambda: s1.stage01_data(force=args.force)),
        2: ("VCF quality control", s1.stage02_qc),
        3: ("coding restriction", s1.stage03_coding),
        4: ("annotation", lambda: s4.stage04_annotate(force=args.force)),
        5: ("recessive scan", s4.stage05_scan),
        6: ("phenotype ranking + permutations", s6.stage06_phenotype),
        7: ("candidate verification", s7.stage07_verify),
        8: ("[supp] locus review", s8.stage08_locus),
        9: ("[supp] NMD position", s9.stage09_nmd),
        10: ("[supp] gnomAD lookup", s9.stage10_gnomad),
        11: ("[supp] allele-fraction screen", s11.stage11_baf),
        12: ("submission CSV", s12.stage12_submission),
        13: ("[test] weight and threshold sensitivity", s13.stage13_robustness),
        14: ("[test] in-silico recovery benchmark", s13.stage14_benchmark),
    }
    timings = {}
    for n in parse_stages(args.stages):
        name, fn = registry[n]
        t0 = time.time()
        print(f"\n=== Stage {n}: {name} ===", flush=True)
        try:
            res = fn()
        except Exception as e:
            print(f"Stage {n} FAILED: {type(e).__name__}: {e}")
            if n <= 6 or n == 12:
                return 1
            continue
        print(json.dumps(res, indent=2, default=str)[:2500])
        timings[n] = {"name": name, "seconds": round(time.time() - t0, 1)}
        print(f"(stage {n} took {time.time() - t0:.1f} s)")
    if timings:
        prev = {}
        tf = cfg.path("timings")
        if tf.exists():
            prev = json.load(open(tf))
        prev.update({str(k): v for k, v in timings.items()})
        with open(tf, "w") as f:
            json.dump(prev, f, indent=2)
    with open(cfg.workdir() / "environment.json", "w") as f:
        json.dump(versions(), f, indent=2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
